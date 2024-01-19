#' getInVitroCompartments
#'
#' @description
#' Generates a list of values describing an in vitro compartment
#'
#' @param cCells cells concentration (in million/ml)
#' @param cSerum serum concentration (in ???)
#' @param microplateWells number of wells in the microplate
#' @param volMedium volume of medium in the well (in mL)
#'
#' @return a list of values representing
#' @export
#'
#' @examples
getInVitroCompartment <- function(cCells, cSerum, microplateType, volMedium) {
  # cCells do have to be in nmillion/ml?
  # the units of the concentrations and surface of plastic are related to how the partitions coefficient were derived

  # calculate the concentration of lipids in the system
  # Based on DOI 10.1002/jps.23602
  # This approximations lead to similar values to what FFischer report for HepG2
  Fapl <- 0.0001 * cCells
  Fnpl <- 5.26 * Fapl # neutral phospholipids in liver is 5.26 times the amount of acidic phospholipids
  Fpl <- Fnpl + Fapl
  Fnl <- 0.000652 * cCells # this includes all neutral lipids ( storage and neutral phospholipids)
  Fw <- 1 - Fnl - Fpl
  cLipCells <- Fnl

  cLipMedium <- 0.00018
  cLipSerum <- cSerum * 0.00157
  # cLipMediumt value is from FFischer 2017 average medium lipid composition
  # the value multiplying with the cSerum is from FFischer 2017 average FBS lipid composition

  cLip <- cLipSerum + cLipCells + cLipMedium

  # calculate the concentration of protein in the system
  cProMedium <- 0.00105
  cProSerum <- cSerum * 0.053 # F Fischer measured 53 mL/L protein content in FBS
  cProCells <- 2.87E-4
  cPro <- cProMedium + cProCells + cProSerum

  # calculate the plastic in the system
  if (microplateType == 96) {
    diam_mm <- 6.6
    volWell_cm3 <- 0.392
  }

  areaGrowth_mm2 <- pi * (diam_mm / 2)^2
  volWell_mm3 <- volWell_cm3 * 1000
  volMedium_mm3 <- volMedium * 1000
  volMedium_L <- volMedium / 1000
  heighMedium_mm3 <- volMedium_mm3 / areaGrowth_mm2
  surfAreaP_mm2 <- 2 * pi * (diam_mm / 2) * heighMedium_mm3
  surfAreaP_m2 <- surfAreaP_mm2 / 1E6
  volAir_L <- (volWell_mm3 - volMedium_mm3) / 1E6
  saPlasticVolMedium <- surfAreaP_m2 / volMedium_L


  #TODO compute all values below before assigning into list

  inVitroCompartment <-
    list(
      # cCellLipN = cCellLipN,
      # cCellLipNPL = cCellLipNPL,
      # cCellLipAPL = cCellLipAPL,
      # cCellPro = cCellPro,
      # cFBSLipN = cFBSLipN,
      # cFBSLipNPL = cFBSLipNPL,
      # cFBSLipAPL = cFBSLipAPL,
      # cFBSPro = cFBSPro,
      # SaPlasticVolMedium = saPlasticVolMedium,
      volAir_L = volAir_L
    )

  return(inVitroCompartment)
}
