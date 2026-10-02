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
#' `method = "pplfer"`. Descriptors can be obtained, for example, from the fup
#' calculator (<https://drumap.nibiohn.go.jp/fup/>).
#'
#' The PP-LFER method uses the polyparameter LFERs at 37 C for membrane
#' lipid-water (Endo et al. 2011), bovine serum albumin-water (Endo and Goss
#' 2011) and muscle protein-water (Endo et al. 2012) partitioning, for the
#' membrane lipids, albumin and globulins respectively. Each system has two
#' published equations (one using L, one using E); both are calculated and the
#' resulting K values are averaged (arithmetic mean in linear space). The
#' protein equations give K per volume of protein (L/L), which is converted to
#' L/kg using a protein density of 1.35 g/mL. The equations are for the
#' neutral species, ionization is not taken into account.
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
  # membrane lipid-water, Endo et al. (2011)
  membrane_lipid = list(
    c(e = 0, s = -0.93, a = -0.18, b = -3.75, v = 1.73, l = 0.49, c = 0.53),
    c(e = 0.74, s = -0.72, a = 0.11, b = -3.63, v = 3.3, l = 0, c = 0.29)
  ),
  # bovine serum albumin-water, Endo and Goss (2011)
  albumin = list(
    c(e = 0, s = -0.46, a = 0.2, b = -3.18, v = 1.84, l = 0.28, c = 0.48),
    c(e = 0.36, s = -0.26, a = 0.37, b = -3.23, v = 2.82, l = 0, c = 0.27)
  ),
  # muscle protein-water, Endo et al. (2012)
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
