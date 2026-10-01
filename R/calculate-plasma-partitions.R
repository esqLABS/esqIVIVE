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
#' @param ionization Ionization class of up to two ionizable groups, as a
#'   vector of length 2 with `"acid"`, `"base"` or `"neutral"`, for example
#'   `c("acid", "neutral")`.
#' @param pka pKa values of the two ionizable groups, a vector of length 2.
#' @param lfer_e,lfer_b,lfer_a,lfer_s,lfer_v Abraham solute descriptors E, B,
#'   A, S and V (the McGowan volume). Needed when `method = "pplfer"`.
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
#' How ionization is taken into account by the PP-LFER method is still under
#' review.
#'
#' @examples
#' calculate_plasma_partitions(
#'   method = "logp",
#'   lipophilicity = 2,
#'   ionization = c("acid", "neutral"),
#'   pka = c(3, 0)
#' )
#'
#' calculate_plasma_partitions(
#'   method = "pplfer",
#'   lipophilicity = 2,
#'   ionization = c("acid", "neutral"),
#'   pka = c(3, 0),
#'   lfer_e = 1,
#'   lfer_b = 0,
#'   lfer_a = 1.5,
#'   lfer_s = 0.8,
#'   lfer_v = 2
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
  verbose = FALSE
) {
  method <- rlang::arg_match(method, c("logp", "pplfer"))
  if (method == "pplfer") {
    descriptors <- list(
      lfer_e = lfer_e,
      lfer_b = lfer_b,
      lfer_a = lfer_a,
      lfer_s = lfer_s,
      lfer_v = lfer_v
    )
    missing <- names(descriptors)[vapply(descriptors, is.null, logical(1))]
    if (length(missing) > 0) {
      cli::cli_abort("{.val pplfer} needs {.arg {missing}}.")
    }
  }

  X <- .calculate_ionization_factors(ionization, pka)[["ion_factor_plasma"]] #Interstitial tissue

  if (method == "logp") {
    logD <- lipophilicity * 1 / (1 + X)

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
  } else {
    #Add LFER_a

    LFER_Ei <- 0.15 + lfer_e

    LFER_Vi <- -0.0215 + lfer_v

    LFER_Bi <- 2.15 - 0.204 * lfer_s + 1.217 * lfer_b + 0.314 * lfer_v

    LFER_Si <- 1.224 + 0.908 * lfer_e + 0.827 * lfer_s + 0.453 * lfer_v

    LFER_Ai <- -0.208 - 0.058 * lfer_s + 0.0354 * lfer_a + 0.076 * lfer_v

    LFER_J <- 1.793 + 0.267 * lfer_e - 0.195 * lfer_s + 0.35 * lfer_v

    #check possibly appli limit, range chemicals...
    kmemlip_LL <- 10^(0.29 +
      0.74 * lfer_e -
      0.72 * lfer_s -
      3.63 * lfer_b +
      3.3 * lfer_v)
    #equation for ions from https://pubs.acs.org/doi/10.1021/acs.est.5b06176
    kbsa_Lkg <- 0.85 +
      0.63 * LFER_Ei -
      0.63 * LFER_Si -
      0.05 * LFER_Ai +
      2.08 * LFER_Bi +
      2.06 * LFER_Vi +
      3.13 * LFER_J
    kalb_Lkg <- kbsa_Lkg
    kmus_Lkg <- -0.24 +
      0.68 * lfer_e -
      0.76 * lfer_s -
      2.29 * lfer_b +
      2.51 * lfer_v
    kglob_Lkg <- kmus_Lkg
  }
  result <- list(
    partition_albumin = kalb_Lkg,
    partition_globulin = kglob_Lkg,
    partition_membrane_lipids = kmemlip_LL
  )

  if (verbose) {
    .print_ivive_result(
      "calculate_plasma_partitions",
      inputs = list(
        method = method,
        lipophilicity = lipophilicity,
        ionization = ionization,
        pka = pka
      ),
      result = result
    )
  }

  result
}
