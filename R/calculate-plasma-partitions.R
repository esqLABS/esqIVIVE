#' Calculate the partition coefficients to plasma components
#'
#' @description
#' Predicts the partition coefficients of a compound to the albumin, the
#' globulins and the membrane lipids (such as those of lipoproteins) of
#' plasma, from its lipophilicity or from its PP-LFER descriptors. Pass the
#' results to [calculate_fu_plasma()] to predict the fraction unbound in
#' plasma.
#'
#' @param method Prediction method: `"logp"` (regressions on lipophilicity)
#'   or `"pplfer"` (poly-parameter linear free energy relationships).
#' @param lipophilicity Lipophilicity of the compound as logP (log units).
#'   Only used when `method = "logp"`.
#' @param ionization Ionization class of up to two ionizable groups, as a
#'   vector of length 2 with `"acid"`, `"base"` or `"neutral"`, for example
#'   `c("acid", "neutral")`. Only used when `method = "logp"`.
#' @param pka pKa values of the two ionizable groups, a vector of length 2.
#'   Only used when `method = "logp"`.
#' @param lfer_e,lfer_b,lfer_a,lfer_s,lfer_v,lfer_l Abraham solute descriptors
#'   E (excess molar refraction), B (hydrogen-bond basicity), A
#'   (hydrogen-bond acidity), S (dipolarity/polarizability), V (McGowan
#'   volume) and L (log of the gas-hexadecane partition coefficient at 25 C).
#'   Needed when `method = "pplfer"`.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @return A named list: `partition_albumin` (L/kg), `partition_globulin`
#'   (L/kg) and `partition_membrane_lipids` (L/L).
#'
#' @details
#' For neutral compounds with logP above 4, use `method = "logp"`. For acidic
#' phenols, carboxylic acids, pyridines and amines you can use
#' `method = "pplfer"`. Descriptors can be obtained, for example, from the
#' UFZ-LSER database (<https://web.app.ufz.de/compbc/lserd/public/start/>) or
#' from the fup calculator (<https://drumap.nibiohn.go.jp/fup/>).
#'
#' The PP-LFER method uses the polyparameter LFERs at 37 C for membrane
#' lipid-water (Endo et al. 2011), bovine serum albumin-water (Endo and Goss
#' 2011) and muscle protein-water (Endo et al. 2012) partitioning, for the
#' membrane lipids, albumin and globulins respectively. The system parameters
#' and the algorithms (two equations per system, one using L and one using E)
#' are those of the UFZ-LSER database, which compiles the published equations;
#' the equations that use L are also given in the supporting information of
#' Endo et al. (2013). Both equations are calculated and the resulting K
#' values are averaged (arithmetic mean in linear space). The equations are for
#' the neutral species, ionization is not taken into account.
#'
#' The protein partition coefficients are published per kg of protein (L/kg).
#' The UFZ-LSER database gives them per volume of protein (L/L), which adds
#' log10(1.35) to the intercept: the intercepts of the muscle protein and
#' albumin equations differ by 0.13 from the papers. The result is converted
#' back to L/kg with a protein density of 1.35 g/mL.
#'
#' The muscle protein equation is the one fitted to chicken muscle protein
#' (Endo et al. 2012), which is similar to the one for fish muscle protein. The
#' albumin equations were determined for bovine serum albumin (heat shock
#' fraction, fatty acid free) and fit worse than the others (SD of 0.41 log
#' units, against 0.22 to 0.28 for the muscle protein and membrane lipid
#' equations), because binding to albumin is partly specific. Binding to the
#' globulins is assumed to be that of muscle protein.
#'
#' @references
#' UFZ-LSER database v 4.1.2. Helmholtz Centre for Environmental
#' Research-UFZ, Leipzig, Germany. <https://www.ufz.de/lserd>. Web
#' application: <https://web.app.ufz.de/compbc/lserd/public/start/>
#'
#' Endo S, Escher BI, Goss K-U (2011). Capacities of membrane lipids to
#' accumulate neutral organic chemicals. *Environmental Science & Technology*
#' 45(14):5912-5921. <https://doi.org/10.1021/es200855w>
#'
#' Endo S, Goss K-U (2011). Serum albumin binding of structurally diverse
#' neutral organic compounds: data and models. *Chemical Research in
#' Toxicology* 24(12):2293-2301. <https://doi.org/10.1021/tx200431b>
#'
#' Endo S, Bauerfeind J, Goss K-U (2012). Partitioning of neutral organic
#' compounds to structural proteins. *Environmental Science & Technology*
#' 46(22):12697-12703. <https://doi.org/10.1021/es303379y>
#'
#' Endo S, Brown TN, Goss K-U (2013). General model for estimating partition
#' coefficients to organisms and their tissues using the biological
#' compositions and polyparameter linear free energy relationships.
#' *Environmental Science & Technology* 47(12):6630-6639.
#' <https://doi.org/10.1021/es401772m>
#'
#' @examples
#' calculate_plasma_partitions(
#'   method = "logp",
#'   lipophilicity = 2,
#'   ionization = c("acid", "neutral"),
#'   pka = c(3, 0)
#' )
#'
#' # Abraham descriptors of ketoconazole
#' calculate_plasma_partitions(
#'   method = "pplfer",
#'   lfer_e = 3.14,
#'   lfer_s = 3.68,
#'   lfer_a = 0,
#'   lfer_b = 2.6,
#'   lfer_v = 3.7208,
#'   lfer_l = 19.71
#' )
#' @export
calculate_plasma_partitions <- function(
  method,
  lipophilicity,
  ionization,
  pka,
  lfer_e = NULL,
  lfer_b = NULL,
  lfer_a = NULL,
  lfer_s = NULL,
  lfer_v = NULL,
  lfer_l = NULL,
  verbose = FALSE
) {
  method <- rlang::arg_match(method, c("logp", "pplfer"))

  if (method == "logp") {
    groups <- .order_ionizable_groups(ionization, pka)
    ionization <- groups$ionization
    pka <- groups$pka

    X <- .calculate_ionization_factors(ionization, pka)[["ion_factor_plasma"]] #Interstitial tissue

    logD <- log10(10^lipophilicity * 1 / (1 + X))

    if (pka[1] != 0) {
      kmemlip_LL <- 10^logD
    } else {
      #Yu et al  regression
      kmemlip_LL <- 10^(1.294 + 0.304 * lipophilicity)
    }

    #for albumin we are not correcting for ionization since acid molecules also bind albumin

    kalb_Lkg <- 0.163 + 0.0221 * kmemlip_LL #Schmitt equation

    kglob_Lkg_1 <- 0.163 + 0.0221 * kmemlip_LL #Schmitt equation for general tissue protein

    kglob_Lkg_2 <- 10^(0.37 * logD - 0.29) #based on the eq used in the VCBA

    kglob_Lkg <- mean(c(kglob_Lkg_1, kglob_Lkg_2))

    inputs <- list(
      method = method,
      lipophilicity = lipophilicity,
      ionization = ionization,
      pka = pka
    )
  } else {
    descriptors <- list(
      lfer_e = lfer_e,
      lfer_b = lfer_b,
      lfer_a = lfer_a,
      lfer_s = lfer_s,
      lfer_v = lfer_v,
      lfer_l = lfer_l
    )
    missing <- names(descriptors)[vapply(descriptors, is.null, logical(1))]
    if (length(missing) > 0) {
      cli::cli_abort("{.val pplfer} needs {.arg {missing}}.")
    }

    kmemlip_LL <- .calculate_pplfer_partition(
      .pplfer_coefficients$membrane_lipid,
      lfer_e,
      lfer_s,
      lfer_a,
      lfer_b,
      lfer_v,
      lfer_l
    )

    # the protein equations give K per volume of protein (L/L), convert to per kg
    protein_density_kgL <- 1.35 # g/mL
    kalb_Lkg <- .calculate_pplfer_partition(
      .pplfer_coefficients$albumin,
      lfer_e,
      lfer_s,
      lfer_a,
      lfer_b,
      lfer_v,
      lfer_l
    ) /
      protein_density_kgL
    kglob_Lkg <- .calculate_pplfer_partition(
      .pplfer_coefficients$muscle_protein,
      lfer_e,
      lfer_s,
      lfer_a,
      lfer_b,
      lfer_v,
      lfer_l
    ) /
      protein_density_kgL

    inputs <- c(list(method = method), descriptors)
  }
  result <- list(
    partition_albumin = kalb_Lkg,
    partition_globulin = kglob_Lkg,
    partition_membrane_lipids = kmemlip_LL
  )

  if (verbose) {
    .print_ivive_result("calculate_plasma_partitions", inputs = inputs, result = result)
  }

  result
}

# PP-LFER coefficients, logK = c + e*E + s*S + a*A + b*B + v*V + l*L, at 37 C.
# Each system has two published equations: one that uses L (and no E) and one
# that uses E (and no L). Terms missing in the source tables are 0.
.pplfer_coefficients <- list(
  # membrane lipid-water, Endo et al. (2011). The L equation matches Endo et
  # al. (2013) SI Table S1, the E equation is from UFZ-LSER
  membrane_lipid = list(
    c(e = 0, s = -0.93, a = -0.18, b = -3.75, v = 1.73, l = 0.49, c = 0.53),
    c(e = 0.74, s = -0.72, a = 0.11, b = -3.63, v = 3.3, l = 0, c = 0.29)
  ),
  # bovine serum albumin-water, Endo and Goss (2011), Table 2. The intercepts
  # are per volume of protein (UFZ-LSER): 0.13 = log10(1.35) above the paper,
  # which gives 0.35 (L equation) and 0.14 (E equation) in L/kg
  albumin = list(
    c(e = 0, s = -0.46, a = 0.2, b = -3.18, v = 1.84, l = 0.28, c = 0.48),
    c(e = 0.36, s = -0.26, a = 0.37, b = -3.23, v = 2.82, l = 0, c = 0.27)
  ),
  # muscle protein-water, Endo et al. (2012), chicken muscle protein, eq 6 for
  # the E equation and Endo et al. (2013) SI Table S1 for the L equation. The
  # E intercept is per volume of protein (UFZ-LSER): -0.65 vs -0.79 in the paper
  muscle_protein = list(
    c(e = 0, s = -0.59, a = 0.21, b = -3.17, v = 2.13, l = 0.33, c = -0.94),
    c(e = 0.51, s = -0.51, a = 0.26, b = -2.98, v = 3.01, l = 0, c = -0.65)
  )
)

# Calculates logK with each equation and returns the average K (linear space)
.calculate_pplfer_partition <- function(equations, E, S, A, B, V, L) {
  log_k <- vapply(
    equations,
    function(cf) {
      cf[["c"]] +
        cf[["e"]] * E +
        cf[["s"]] * S +
        cf[["a"]] * A +
        cf[["b"]] * B +
        cf[["v"]] * V +
        cf[["l"]] * L
    },
    numeric(1)
  )
  mean(10^log_k)
}
