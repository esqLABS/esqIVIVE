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

# Give a warning for the scaling factors that are NA in the table (no value
# available for this species, tissue and system). The scaling factor, and so the
# result, is then NA.
.warn_unsupported_scaling_factors <- function(
  scaling_factors,
  columns,
  species,
  tissue,
  call = rlang::caller_env()
) {
  descriptions <- c(
    fcell = "the fraction of cells in the tissue",
    weightorgankgBW = "the organ weight per kg body weight",
    MicProtGO = "the microsomal protein per gram tissue",
    CytosProtGO = "the cytosolic protein per gram tissue",
    CellsGO = "the cells per gram tissue"
  )
  for (column in columns) {
    if (is.na(scaling_factors[[column]])) {
      cli::cli_warn(
        c(
          "The scaling factor {.field {column}} ({descriptions[[column]]}) is
           not supported for species {.val {species}} and tissue
           {.val {tissue}}: it is {.code NA} in the scaling factor table.",
          "i" = "The result is {.code NA}."
        ),
        call = call
      )
    }
  }
  invisible(NULL)
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
