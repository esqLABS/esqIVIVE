#' Scale an in vitro clearance to in vivo
#'
#' @description
#' Converts an in vitro half-life, depletion rate constant or intrinsic
#' clearance into the specific clearance (1/min) that PK-Sim expects. The
#' value is corrected for binding in the incubation and scaled to the whole
#' tissue with the physiological scaling factors of the species.
#'
#' @param value_type Type of `value`:
#'   * `"half_life"`: in vitro half-life of substrate depletion;
#'   * `"rate_constant"`: depletion rate constant, for example from
#'     [fit_depletion_curve()];
#'   * `"intrinsic_clearance"`: in vitro intrinsic clearance.
#' @param value Measured value, in `unit`.
#' @param unit Unit of `value`. The allowed units depend on `value_type`:
#'   * `"half_life"`: `"minutes"`, `"hours"`, `"seconds"`.
#'   * `"rate_constant"`: `"/minutes"`, `"/hours"`, `"/seconds"`.
#'   * `"intrinsic_clearance"` per million cells, per cell or per mg protein:
#'     a volume (`mL`, `uL`, `L`), a time (`minutes`, `hours`, `seconds`) and
#'     the amount, for example `"mL/minutes/millioncells"`, `"uL/hours/cell"`
#'     or `"L/minutes/mg protein"`. `L` is not available per cell.
#'   * `"intrinsic_clearance"` per incubation: `"mL/minutes"`, `"uL/minutes"`,
#'     `"mL/seconds"`, `"uL/seconds"`, `"mL/hours"`, `"uL/hours"`.
#'   * `"intrinsic_clearance"` per kg body weight: `"mL/minutes/kg"`,
#'     `"uL/minutes/kg"`, `"mL/hours/kg"`, `"uL/hours/kg"`.
#' @param system Incubation system: `"microsomes"`, `"cells"` (for example
#'   hepatocytes) or `"cytosol"` (cytosolic fraction). Microsomes are scaled
#'   with the microsomal protein per gram tissue, cells with the cells per gram
#'   tissue and the cytosol with the cytosolic protein per gram tissue.
#' @param fu_in_vitro Fraction unbound in the incubation, for example from
#'   [calculate_fu_in_vitro()]. Defaults to 1 (no binding).
#' @param concentration_microsomes Microsomal protein concentration (mg/mL).
#'   Needed for microsomes when `value_type` is `"half_life"` or
#'   `"rate_constant"`, or when `unit` is per incubation.
#' @param concentration_cells Cell concentration (million cells/mL). Needed
#'   for cells in the same cases.
#' @param concentration_cytosol Cytosolic protein concentration (mg/mL).
#'   Needed for the cytosol in the same cases.
#' @param volume_medium Volume of medium in the incubation (mL), used for the
#'   per incubation units. Defaults to 1.
#' @param empirical_correction If `TRUE`, apply the empirical correction
#'   factors of Wood et al. (2017), which correct the tendency of in vitro data
#'   to overpredict slow and underpredict fast clearances. Available for human
#'   and rat, with microsomes or cells.
#' @param tissue Tissue whose scaling factors are used. Defaults to `"liver"`.
#' @param species Species whose scaling factors are used: `"human"`, `"rat"`
#'   or `"dog"`. Defaults to `"human"`.
#' @param relative_expression_factor Relative expression or activity factor of
#'   the enzyme in vivo compared with the incubation. Defaults to 1.
#' @param verbose If `TRUE`, print the inputs.
#'
#' @return The specific clearance (1/min), a single number.
#' @export
#' @examples
#' # cells, for example hepatocytes
#' ivive_clearance(
#'   value_type = "intrinsic_clearance",
#'   value = 18.27,
#'   unit = "mL/minutes/millioncells",
#'   system = "cells",
#'   fu_in_vitro = 0.5,
#'   concentration_cells = 0.5
#' )
#'
#' # microsomes, from the in vitro half-life
#' ivive_clearance(
#'   value_type = "half_life",
#'   value = 3.9,
#'   unit = "minutes",
#'   system = "microsomes",
#'   fu_in_vitro = 0.4,
#'   concentration_microsomes = 1
#' )
#'
#' # cytosolic fraction
#' ivive_clearance(
#'   value_type = "intrinsic_clearance",
#'   value = 0.05,
#'   unit = "mL/minutes/mg protein",
#'   system = "cytosol",
#'   fu_in_vitro = 0.8,
#'   concentration_cytosol = 1
#' )
ivive_clearance <- function(
  value_type,
  value,
  unit,
  system,
  fu_in_vitro = 1,
  concentration_microsomes = NULL,
  concentration_cells = NULL,
  concentration_cytosol = NULL,
  volume_medium = 1,
  empirical_correction = FALSE,
  tissue = "liver",
  species = "human",
  relative_expression_factor = 1,
  verbose = FALSE
) {
  # check if the arguments are valid
  value_type <- rlang::arg_match(
    value_type,
    c("half_life", "rate_constant", "intrinsic_clearance")
  )
  system <- rlang::arg_match(system, c("microsomes", "cells", "cytosol"))
  unit <- rlang::arg_match(unit, .clearance_units[[value_type]])
  if (!rlang::is_bool(empirical_correction)) {
    cli::cli_abort(
      "{.arg empirical_correction} must be {.code TRUE} or {.code FALSE}, not
       {.val {empirical_correction}}."
    )
  }
  .check_fu_in_vitro_value(fu_in_vitro)

  scaling_factors <- .get_scaling_factors(species, tissue)
  fintcell <- scaling_factors[["fcell"]]
  organkgBW <- scaling_factors[["weightorgankgBW"]]

  #chose the system specific scaling factors
  if (system == "microsomes") {
    nLiver <- scaling_factors[["MicProtGO"]] # mg protein/g liver
    cInvitro <- concentration_microsomes #mg/mL
    concentration <- "concentration_microsomes"
  } else if (system == "cytosol") {
    nLiver <- scaling_factors[["CytosProtGO"]] # mg protein/g liver
    cInvitro <- concentration_cytosol #mg/mL
    concentration <- "concentration_cytosol"
  } else {
    nLiver <- scaling_factors[["CellsGO"]]
    cInvitro <- concentration_cells # million cells/mL assay
    concentration <- "concentration_cells"
  }
  needs_concentration <- value_type != "intrinsic_clearance" ||
    unit %in% names(.clearance_per_incubation)
  if (needs_concentration && is.null(cInvitro)) {
    cli::cli_abort(
      "{.arg {concentration}} is needed for {system} when {.arg value_type} is
       {.val {value_type}} and {.arg unit} is {.val {unit}}."
    )
  }

  #combination scale factors
  SF <- nLiver / fintcell * relative_expression_factor / fu_in_vitro
  SF2 <- 1 / fintcell * relative_expression_factor / fu_in_vitro

  #Derive the in vitro clearance value---------------------
  ClspePermin <- if (value_type == "half_life") {
    multFactorHalf <- c(minutes = 1, hours = 1 / 60, seconds = 60)[[unit]]
    kcat_min <- 0.693 / value * multFactorHalf
    kcat_min / cInvitro * SF
  } else if (value_type == "rate_constant") {
    value * .clearance_rate_constants[[unit]] / cInvitro * SF
  } else if (unit %in% names(.clearance_per_amount)) {
    value * .clearance_per_amount[[unit]] * SF
  } else if (unit %in% names(.clearance_per_incubation)) {
    value *
      .clearance_per_incubation[[unit]] /
      (cInvitro * volume_medium) *
      SF
  } else {
    #multiplying by the g liver weight per kilo BW
    value * .clearance_per_kg[[unit]] / organkgBW * SF2
  }

  if (empirical_correction) {
    ClspePermin <- .apply_wood_correction(
      ClspePermin,
      organkgBW = organkgBW,
      system = system,
      species = species
    )
  }

  if (verbose) {
    .print_ivive_result(
      "ivive_clearance",
      inputs = list(
        value_type = value_type,
        value = value,
        unit = unit,
        system = system,
        fu_in_vitro = fu_in_vitro,
        empirical_correction = empirical_correction,
        tissue = tissue,
        species = species
      )
    )
  }

  ClspePermin
}

# Unit conversion factors of ivive_clearance(), before the scaling factors.
.clearance_rate_constants <- c(
  "/minutes" = 1,
  "/hours" = 1 / 60,
  "/seconds" = 60
)

.clearance_per_amount <- c(
  "mL/minutes/millioncells" = 1,
  "uL/minutes/millioncells" = 1 / 1000,
  "L/minutes/millioncells" = 1000,
  "mL/hours/millioncells" = 1 / 60,
  "uL/hours/millioncells" = 1 / 60000,
  "L/hours/millioncells" = 1000 / 60,
  "mL/seconds/millioncells" = 60,
  "uL/seconds/millioncells" = 60 / 1000,
  "L/seconds/millioncells" = 60000,
  "mL/minutes/cell" = 1E6,
  "uL/minutes/cell" = 1000,
  "mL/hours/cell" = 1E6 / 60,
  "uL/hours/cell" = 1E6 / 60000,
  "mL/seconds/cell" = 60 * 1E6,
  "uL/seconds/cell" = 60 * 1E6 / 1000,
  "mL/minutes/mg protein" = 1,
  "uL/minutes/mg protein" = 1 / 1000,
  "L/minutes/mg protein" = 1000,
  "mL/hours/mg protein" = 1 / 60,
  "uL/hours/mg protein" = 1 / 60000,
  "L/hours/mg protein" = 1000 / 60,
  "mL/seconds/mg protein" = 60,
  "uL/seconds/mg protein" = 60 / 1000,
  "L/seconds/mg protein" = 60000
)

.clearance_per_incubation <- c(
  "mL/minutes" = 1,
  "uL/minutes" = 1 / 1000,
  "mL/seconds" = 60,
  "uL/seconds" = 60 / 1000,
  "mL/hours" = 1 / 60,
  "uL/hours" = 1 / 60000
)

.clearance_per_kg <- c(
  "mL/minutes/kg" = 1,
  "uL/minutes/kg" = 1 / 1000,
  "mL/hours/kg" = 1 / 60,
  "uL/hours/kg" = 1 / 60000
)

.clearance_units <- list(
  half_life = c("minutes", "hours", "seconds"),
  rate_constant = names(.clearance_rate_constants),
  intrinsic_clearance = c(
    names(.clearance_per_amount),
    names(.clearance_per_incubation),
    names(.clearance_per_kg)
  )
)

#Add scalars from Wood et al 2017-https://doi.org/10.1124/dmd.117.077040.
# this scalar are empitical and they are based on in vitro data tending to overestimate slow clearance
# and underpredict fast clearance
.apply_wood_correction <- function(
  ClspePermin,
  organkgBW,
  system,
  species,
  call = rlang::caller_env()
) {
  wood_sf <- list()
  wood_sf[["human"]] <- data.frame(
    Cl_ranges = c("<10", "10-100", "100-1000", "1000-10000", ">10000"),
    microsomes = c(0.7, 1.8, 4.6, 7.5, 58),
    cells = c(0.61, 3.9, 7.1, 22, 1200)
  )
  wood_sf[["rat"]] <- data.frame(
    Cl_ranges = c("<10", "10-100", "100-1000", "1000-10000", ">10000"),
    microsomes = c(0.086, 0.83, 1.7, 2.5, 230),
    cells = c(0.13, 1.6, 3.2, 7.2, 180)
  )
  if (!species %in% names(wood_sf)) {
    cli::cli_abort(
      "The empirical correction is only available for {names(wood_sf)}, not
       {.val {species}}.",
      call = call
    )
  }
  if (!system %in% colnames(wood_sf[[species]])) {
    cli::cli_abort(
      "The empirical correction is only available for microsomes and cells,
       not {.val {system}}.",
      call = call
    )
  }
  wood_table <- wood_sf[[species]]

  # conditional to pick the right wood scalar
  scaled <- ClspePermin * organkgBW
  band <- if (scaled < 10) {
    "<10"
  } else if (scaled < 100) {
    "10-100"
  } else if (scaled < 1000) {
    "100-1000"
  } else if (scaled < 10000) {
    "1000-10000"
  } else {
    ">10000"
  }
  ClspePermin * wood_table[wood_table$Cl_ranges == band, system]
}
