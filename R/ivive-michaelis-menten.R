#' Scale Michaelis-Menten parameters to in vivo
#'
#' @description
#' Scales an in vitro Vmax to the whole tissue with the physiological scaling
#' factors of the species, and corrects Km for binding in the incubation.
#'
#' @param system Incubation system, `"microsomes"` or `"cells"` (for
#'   example hepatocytes).
#' @param vmax In vitro Vmax (umol/min/million cells for cells,
#'   umol/min/mg protein for microsomes), for example from
#'   [fit_michaelis_menten_curve()].
#' @param km In vitro Km (uM).
#' @param fu_in_vitro Fraction unbound in the incubation, for example from
#'   [calculate_fu_in_vitro()]. Defaults to 1 (no binding).
#' @param tissue Tissue whose scaling factors are used. Defaults to `"liver"`.
#' @param species Species whose scaling factors are used: `"human"`, `"rat"`
#'   or `"dog"`. Defaults to `"human"`.
#' @param relative_expression_factor Relative expression or activity factor of
#'   the enzyme in vivo compared with the incubation. Defaults to 1. To use it,
#'   set the reference concentration of the enzyme in PK-Sim to 1 uM.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @returns A named list:
#'   * `vmax`: in vivo Vmax (umol/min/L of tissue);
#'   * `km_unbound`: unbound Km (uM).
#' @export
#'
#' @examples
#' ivive_michaelis_menten(system = "cells", vmax = 2, km = 1)
#'
#' ivive_michaelis_menten(
#'   system = "microsomes",
#'   vmax = 2,
#'   km = 1,
#'   fu_in_vitro = 0.2
#' )
ivive_michaelis_menten <- function(
  system,
  vmax,
  km,
  fu_in_vitro = 1,
  tissue = "liver",
  species = "human",
  relative_expression_factor = 1,
  verbose = FALSE
) {
  # check if the arguments are valid
  system <- rlang::arg_match(system, c("microsomes", "cells"))
  .check_fu_in_vitro_value(fu_in_vitro)

  #Correct Km for fraction unbound
  Km_unb_uM <- km * fu_in_vitro

  #Calculate in vivo Vmax--------------------------------------------------------
  scaling_factors <- .get_scaling_factors(species, tissue)
  fintcell <- scaling_factors[["fcell"]]

  #chose the system specific scaling factors
  if (system == "microsomes") {
    scfactor <- scaling_factors[["MicProtGO"]] # mg protein/g liver
  } else {
    scfactor <- scaling_factors[["CellsGO"]]
  }
  dens <- 1000 #g/L
  vmax_umol_minL <- vmax *
    scfactor *
    relative_expression_factor /
    fintcell *
    dens

  result <- list(vmax = vmax_umol_minL, km_unbound = Km_unb_uM)

  if (verbose) {
    .print_ivive_result(
      "ivive_michaelis_menten",
      inputs = list(
        system = system,
        vmax = vmax,
        km = km,
        fu_in_vitro = fu_in_vitro,
        tissue = tissue,
        species = species,
        relative_expression_factor = relative_expression_factor
      ),
      result = result
    )
  }

  result
}
