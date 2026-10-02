#' Calculate the intestinal transcellular permeability
#'
#' @description
#' Converts a measured permeability into the intestinal transcellular
#' permeability (Pint) used by PK-Sim, with an empirical regression calibrated
#' on compounds of high solubility.
#'
#' @param method Type of the measured permeability: `"caco2"` for the
#'   apparent permeability in Caco-2 cells (Papp), or `"peff"` for the human
#'   effective intestinal permeability (Peff), measured in vivo or predicted.
#' @param permeability Measured permeability (cm/s).
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @returns The intestinal transcellular permeability (cm/s), a single number.
#' @export
#'
#' @examples
#' calculate_pint(method = "caco2", permeability = 2.3e-6)
#'
#' calculate_pint(method = "peff", permeability = 2.3e-6)
calculate_pint <- function(method, permeability, verbose = FALSE) {
  method <- rlang::arg_match(method, c("caco2", "peff"))

  #based on calibration with high solubility
  pint_cms <- switch(
    method,
    caco2 = 0.0001 * 10^((0.4428 * (log10(permeability * 1000000) - 3.0941))),
    peff = 0.0001 * 10^((0.864 * (log10(permeability * 1000) - 3.1029)))
  )

  if (verbose) {
    .print_ivive_result(
      "calculate_pint",
      inputs = list(method = method, permeability = permeability),
      result = pint_cms
    )
  }

  pint_cms
}
