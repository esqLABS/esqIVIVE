#' getInVitroCompartments
#'
#' @description
#' Generates a list of values describing an in vitro compartment
#' This function is used inside the Fraction unbound function
#' @param typeSystem if system is hepatocytes or microsomes
#' @param cCells_Mml cells concentration (in million/ml)
#' @param cMicro concentration of microsome protein mg/mL
#' @param FBS_fraction fraction of serum concentration, values can only go from 0-1
#' @param microplateWells number of wells in the microplate
#' @param volMedium volume of medium in the well (in mL)
#'
#' @return a list of values representing the  different in vitro compartments, concentrations are given as fraction of volume
#' @export
#'
#' @examples 
#' getInVitroCompartment("hepatocytes", FBS_fraction=0.05, microplateType = 96, volMedium_mL = 0.15, cCells_Mml = 0.1)
#' getInVitroCompartment("microsomes", FBS_fraction=0, microplateType = 24, volMedium_mL = 0.5, cMicro_mgml = 1)
#'
getInVitroCompartment <- function(
  typeSystem,
  FBS_fraction,
  microplateType,
  volMedium_mL,
  cCells_Mml = NULL,
  cMicro_mgml = NULL
) {
  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes", "microsomes"))

  # Plastic and headspace compartments------------------------------------------
  # calculate the plastic in the system
  if (microplateType == 96) {
    diam_mm <- 6.6
    volWell_cm3 <- 0.392
  } else if (microplateType == 48) {
    diam_mm <- 11
    volWell_cm3 <- 1.62
  } else if (microplateType == 24) {
    diam_mm <- 15.55
    volWell_cm3 <- 3.47
  } else if (microplateType == 12) {
    diam_mm <- 22
    volWell_cm3 <- 6.9
  }

  # the units of the concentrations and surface of plastic are related
  # to how the partitions coefficient were derived
  areaGrowth_mm2 <- pi * (diam_mm / 2)^2
  volWell_mm3 <- volWell_cm3 * 1000
  volMedium_mm3 <- volMedium_mL * 1000
  volMedium_L <- volMedium_mL / 1000
  heighMedium_mm3 <- volMedium_mm3 / areaGrowth_mm2
  surfAreaP_mm2 <- 2 * pi * (diam_mm / 2) * heighMedium_mm3
  surfAreaP_m2 <- surfAreaP_mm2 / 1E6
  volAir_L <- (volWell_mm3 - volMedium_mm3) / 1E6
  saPlasticVolMedium_m2L <- surfAreaP_m2 / volMedium_L

  # Protein and lipid compartments----------------------------------------------
  # density of lipids assumed 0.9 g/mL and of proteins 1.35 g/mL

  if (typeSystem == "hepatocytes") {
    # see report on input parameters for refernces of values
    # these values are going to be lower than Poulin paper indicates
    cellVol_mLM <- 0.00254 # mL per million cells
    cCellAPL_vvmedium <- 0.0088 * cellVol_mLM * cCells_Mml
    cCellNPL_vvmedium <- 0.0331 * cellVol_mLM * cCells_Mml
    cCellPL_vvmedium <- cCellNPL_vvmedium + cCellAPL_vvmedium
    cCellNL_vvmedium <- 0.0445 * cellVol_mLM * cCells_Mml # this includes all neutral lipids ( storage and neutral phospholipids)
    cCellPro_vvcell <- 0.2
    cCellPro_vvmedium <- cCellPro_vvcell * cellVol_mLM * cCells_Mml
  } else if (typeSystem == "microsomes") {
    #despite Poulin showing rat and human separatly it does not appear there issignificant differences
    cCellPro_vvmedium <- cMicro_mgml / 1000 / 1.35 # mg to g to ml
    cCellPL_mgPLmgprot <- 0.797
    cCellPL_vvmedium <- cCellPL_mgPLmgprot * cMicro_mgml / 1000 / 0.9
    cCellNL_mgPLprot <- 0.235
    cCellNL_vvmedium <- cCellNL_mgPLprot * cMicro_mgml / 1000 / 0.9
    cCellAPL_vPLvNL <- 0.18
    cCellAPL_vvmedium <- cCellAPL_vPLvNL * cCellPL_vvmedium
    cCellNPL_vvmedium <- cCellPL_vvmedium - cCellAPL_vvmedium
    # for microsome system consider assay is performed in glass
    saPlasticVolMedium_m2L <- 0
  } else {}

  cMediumNL_vvmedium <- FBS_fraction * 0.00157
  cMediumNPL_vvmedium <- FBS_fraction * 0.0003
  # the value multiplying with the cSerum is from FFischer 2017 average FBS lipid composition

  # calculate the concentration of protein in the system
  # From average protein content in medium from Fischer paper
  cMediumPro_vvmedium <- FBS_fraction * 0.040

  inVitroCompartment <-
    list(
      cCellNL_vvmedium = cCellNL_vvmedium,
      cCellNPL_vvmedium = cCellNPL_vvmedium,
      cCellAPL_vvmedium = cCellAPL_vvmedium,
      cCellPro_vvmedium = cCellPro_vvmedium,
      cMediumNL_vvmedium = cMediumNL_vvmedium,
      cMediumNPL_vvmedium = cMediumNPL_vvmedium,
      cMediumPro_vvmedium = cMediumPro_vvmedium,
      saPlasticVolMedium_m2L = saPlasticVolMedium_m2L,
      volAir_L = volAir_L
    )

  return(inVitroCompartment)
}
