#' Calculate the compartments of an in vitro incubation
#'
#' @description
#' Describes an incubation of cells (such as hepatocytes) or microsomes in a
#' microplate well: the lipid and protein content of the cells or microsomes
#' and of the serum in the medium, the plastic surface in contact with the
#' medium, and the air above it. [calculate_fu_in_vitro()] uses these values for its partition
#' models. You can also use them to describe a virtual incubation.
#'
#' @param system Incubation system, `"microsomes"` or `"cells"` (for
#'   example hepatocytes).
#' @param fbs_fraction Fraction of fetal bovine serum in the medium (0 to 1).
#' @param microplate_type Number of wells of the microplate: 96, 48, 24 or 12.
#' @param volume_medium Volume of medium in the well (mL).
#' @param concentration_cells Cell concentration (million cells/mL).
#'   Needed when `system = "cells"`.
#' @param concentration_microsomes Microsomal protein concentration (mg/mL).
#'   Needed when `system = "microsomes"`.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @return A named list. Lipid and protein contents are fractions of the
#'   medium volume (L/L).
#'   * `cell_neutral_lipids`, `cell_neutral_phospholipids`,
#'     `cell_acidic_phospholipids`, `cell_proteins`: content of the cells or
#'     microsomes.
#'   * `medium_neutral_lipids`, `medium_neutral_phospholipids`,
#'     `medium_proteins`: content of the serum in the medium.
#'   * `plastic_area_per_volume`: plastic surface in contact with the medium
#'     per medium volume (m2/L). Zero for microsomes, which are assumed to be
#'     incubated in glass.
#'   * `volume_air`: volume of air above the medium in the well (L).
#' @export
#'
#' @examples
#' calculate_in_vitro_compartments(
#'   system = "cells",
#'   fbs_fraction = 0.05,
#'   microplate_type = 96,
#'   volume_medium = 0.15,
#'   concentration_cells = 0.1
#' )
#'
#' calculate_in_vitro_compartments(
#'   system = "microsomes",
#'   fbs_fraction = 0,
#'   microplate_type = 24,
#'   volume_medium = 0.5,
#'   concentration_microsomes = 1
#' )
calculate_in_vitro_compartments <- function(
  system,
  fbs_fraction,
  microplate_type,
  volume_medium,
  concentration_cells = NULL,
  concentration_microsomes = NULL,
  verbose = FALSE
) {
  system <- rlang::arg_match(system, c("microsomes", "cells"))
  if (system == "microsomes" && is.null(concentration_microsomes)) {
    cli::cli_abort("{.arg concentration_microsomes} is needed for microsomes.")
  }
  if (system == "cells" && is.null(concentration_cells)) {
    cli::cli_abort("{.arg concentration_cells} is needed for cells.")
  }

  compartments <- c(
    .calculate_cell_compartments(
      system,
      concentration_cells = concentration_cells,
      concentration_microsomes = concentration_microsomes
    ),
    .calculate_medium_compartments(fbs_fraction),
    .calculate_well_geometry(system, microplate_type, volume_medium)
  )

  if (verbose) {
    .print_ivive_result(
      "calculate_in_vitro_compartments",
      inputs = list(
        system = system,
        fbs_fraction = fbs_fraction,
        microplate_type = microplate_type,
        volume_medium = volume_medium,
        concentration_cells = concentration_cells,
        concentration_microsomes = concentration_microsomes
      ),
      result = compartments
    )
  }

  compartments
}

# Protein and lipid content of the cells or microsomes, as fractions of the
# medium volume. Density of lipids assumed 0.9 g/mL and of proteins 1.35 g/mL.
.calculate_cell_compartments <- function(
  system,
  concentration_cells = NULL,
  concentration_microsomes = NULL
) {
  if (system == "cells") {
    # see report on input parameters for refernces of values
    # these values are going to be lower than Poulin paper indicates
    cellVol_mLM <- 0.00254 # mL per million cells
    cCellAPL_vvmedium <- 0.0088 * cellVol_mLM * concentration_cells
    cCellNPL_vvmedium <- 0.0331 * cellVol_mLM * concentration_cells
    cCellNL_vvmedium <- 0.0445 * cellVol_mLM * concentration_cells # this includes all neutral lipids ( storage and neutral phospholipids)
    cCellPro_vvcell <- 0.2
    cCellPro_vvmedium <- cCellPro_vvcell * cellVol_mLM * concentration_cells
  } else {
    #despite Poulin showing rat and human separatly it does not appear there is significant differences
    cCellPro_vvmedium <- concentration_microsomes / 1000 / 1.35 # mg to g to ml
    cCellPL_mgPLmgprot <- 0.797
    cCellPL_vvmedium <- cCellPL_mgPLmgprot *
      concentration_microsomes /
      1000 /
      0.9
    cCellNL_mgPLprot <- 0.235
    cCellNL_vvmedium <- cCellNL_mgPLprot * concentration_microsomes / 1000 / 0.9
    cCellAPL_vPLvNL <- 0.18
    cCellAPL_vvmedium <- cCellAPL_vPLvNL * cCellPL_vvmedium
    cCellNPL_vvmedium <- cCellPL_vvmedium - cCellAPL_vvmedium
  }

  list(
    cell_neutral_lipids = cCellNL_vvmedium,
    cell_neutral_phospholipids = cCellNPL_vvmedium,
    cell_acidic_phospholipids = cCellAPL_vvmedium,
    cell_proteins = cCellPro_vvmedium
  )
}

# Protein and lipid content of the serum in the medium, as fractions of the
# medium volume.
.calculate_medium_compartments <- function(fbs_fraction) {
  list(
    # the value multiplying with the cSerum is from FFischer 2017 average FBS lipid composition
    medium_neutral_lipids = fbs_fraction * 0.00157,
    medium_neutral_phospholipids = fbs_fraction * 0.0003,
    # From average protein content in medium from Fischer paper
    medium_proteins = fbs_fraction * 0.040
  )
}

# Plastic surface in contact with the medium (m2/L of medium) and air volume
# in the well (L).
.calculate_well_geometry <- function(
  system,
  microplate_type,
  volume_medium,
  call = rlang::caller_env()
) {
  wells <- list(
    "96" = c(diam_mm = 6.6, volWell_cm3 = 0.392),
    "48" = c(diam_mm = 11, volWell_cm3 = 1.62),
    "24" = c(diam_mm = 15.55, volWell_cm3 = 3.47),
    "12" = c(diam_mm = 22, volWell_cm3 = 6.9)
  )
  if (
    length(microplate_type) != 1 ||
      !as.character(microplate_type) %in% names(wells)
  ) {
    cli::cli_abort(
      "{.arg microplate_type} must be one of 96, 48, 24 or 12, not
       {.val {microplate_type}}.",
      call = call
    )
  }
  diam_mm <- wells[[as.character(microplate_type)]][["diam_mm"]]
  volWell_cm3 <- wells[[as.character(microplate_type)]][["volWell_cm3"]]

  # the units of the concentrations and surface of plastic are related
  # to how the partitions coefficient were derived
  areaGrowth_mm2 <- pi * (diam_mm / 2)^2
  volWell_mm3 <- volWell_cm3 * 1000
  volMedium_mm3 <- volume_medium * 1000
  volMedium_L <- volume_medium / 1000
  heighMedium_mm3 <- volMedium_mm3 / areaGrowth_mm2
  surfAreaP_mm2 <- 2 * pi * (diam_mm / 2) * heighMedium_mm3
  surfAreaP_m2 <- surfAreaP_mm2 / 1E6
  volAir_L <- (volWell_mm3 - volMedium_mm3) / 1E6
  saPlasticVolMedium_m2L <- surfAreaP_m2 / volMedium_L

  if (system == "microsomes") {
    # for microsome system consider assay is performed in glass
    saPlasticVolMedium_m2L <- 0
  }

  list(
    plastic_area_per_volume = saPlasticVolMedium_m2L,
    volume_air = volAir_L
  )
}
