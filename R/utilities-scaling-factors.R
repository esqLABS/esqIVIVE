# Physiological scaling factors for a species and tissue, from
# inst/extdata/scaling_factors.csv: one row with fcell, weightorgankgBW,
# MicProtGO, CytosProtGO and CellsGO.
.get_scaling_factors <- function(species, tissue, call = rlang::caller_env()) {
  path <- system.file("extdata", "scaling_factors.csv", package = "ESQivive")
  scaling_factors <- utils::read.csv(path)

  rlang::arg_match(
    species,
    unique(scaling_factors[, "species"]),
    error_arg = "species",
    error_call = call
  )
  rlang::arg_match(
    tissue,
    unique(scaling_factors[, "organ"]),
    error_arg = "tissue",
    error_call = call
  )

  scaling_factors[
    scaling_factors[, "species"] == species &
      scaling_factors[, "organ"] == tissue,
  ]
}

.check_fu_in_vitro_value <- function(fu_in_vitro, call = rlang::caller_env()) {
  if (!is.numeric(fu_in_vitro) || any(fu_in_vitro <= 0 | fu_in_vitro > 1)) {
    cli::cli_abort(
      "{.arg fu_in_vitro} must be above 0 and at most 1, not
       {.val {fu_in_vitro}}.",
      call = call
    )
  }
  invisible(NULL)
}
