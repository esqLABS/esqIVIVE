#' getInVitroCompartments
#'
#' @description
#' Generates a list of values describing an in vitro compartment
#' This function is used inside the Fraction unbound function
#' @param typeSystem if system is hepatocytes or microsomes
#' @param cCells cells concentration (in million/ml)
#' @param cMicro concentration of microsome protein mg/mL
#' @param FBS fraction of serum concentration, values can only go from 0-1
#' @param microplateWells number of wells in the microplate
#' @param volMedium volume of medium in the well (in mL)
#'
#' @return a list of values representing the  different in vitro compartments, concentrations are given as fraction of volume
#' @export a list of concentration of lipids, volume of headspace and surface area of plastic
#'
#' @examples getInVitroCompartment("hepatocytes",0.05,microplateType=96,volMedium=0.15,cCells=0.1)
#' @examples getInVitroCompartment("microsomes",0,microplateType=24,volMedium=0.5,cMicro=1)
#'
getInVitroCompartment <- function(typeSystem,FBS,microplateType, volMedium,cCells=NULL,cMicro=NULL) {

  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes", "microsomes"))

  #Plastic and headspace compartments------------------------------------------
  # calculate the plastic in the system
  if (microplateType == 96) {
    diam_mm <- 6.6
    volWell_cm3 <- 0.392
  } else if (microplateType == 48){
    diam_mm <- 11
    volWell_cm3 <- 1.62
  } else if (microplateType == 24){
    diam_mm <- 15.55
    volWell_cm3 <- 3.47
  }else if (microplateType == 12){
    diam_mm <- 22
    volWell_cm3 <- 6.9
  }

  # the units of the concentrations and surface of plastic are related
  # to how the partitions coefficient were derived
  areaGrowth_mm2 <- pi * (diam_mm / 2)^2
  volWell_mm3 <- volWell_cm3 * 1000
  volMedium_mm3 <- volMedium * 1000
  volMedium_L <- volMedium / 1000
  heighMedium_mm3 <- volMedium_mm3 / areaGrowth_mm2
  surfAreaP_mm2 <- 2 * pi * (diam_mm / 2) * heighMedium_mm3
  surfAreaP_m2 <- surfAreaP_mm2 / 1E6
  volAir_L <- (volWell_mm3 - volMedium_mm3) / 1E6
  saPlasticVolMedium <- surfAreaP_m2 / volMedium_L

  #Protein and lipod compartments----------------------------------------------
  #density of lipids assumed 0.9 g/mL and of proteins 1.35 g/mL

  if (typeSystem=="hepatocytes"){
  # see report on input parameters for refernces of values
  cellVol<-0.00254 # mL per million cells
  cCellAPL<- 0.0088 * cellVol*cCells
  cCellNPL<-0.0331* cellVol*cCells
  cCellPL <- cCellNPL + cCellAPL
  cCellNL <- 0.0445* cellVol*cCells# this includes all neutral lipids ( storage and neutral phospholipids)
  cCellPro <- 0.2*0.00254*cCells

  } else if (typeSystem=="microsomes"){

  cCellPro<- cMicro/1000/1.35 #mg to g to ml
  cCellPL<- 0.797 * cMicro /1000 /0.9
  cCellNL<- 0.235 * cMicro /1000 /0.9 # 0.235 mg lipid/mg protein
  cCellAPL<- 0.18 * cCellPL
  cCellNPL<- cCellPL- cCellAPL
  #for microsome system consider assay is performed in glass
  saPlasticVolMedium=0
  } else {}

  cMediumNL <- FBS * 0.00157
  cMediumNPL <- FBS * 0.0003
  # the value multiplying with the cSerum is from FFischer 2017 average FBS lipid composition

  # calculate the concentration of protein in the system
  # From average protein content in medium from Fischer paper
  cMediumPro <- FBS * 0.040

  inVitroCompartment <-
    list(
      cCellNL = cCellNL,
      cCellNPL = cCellNPL,
      cCellAPL = cCellAPL,
      cCellPro = cCellPro,
      cMediumNL = cMediumNL,
      cMediumNPL = cMediumNPL,
      cMediumPro = cMediumPro,
      saPlasticVolMedium = saPlasticVolMedium,
      volAir_L = volAir_L
    )

  return(inVitroCompartment)
}
