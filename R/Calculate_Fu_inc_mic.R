#' getInVitroFractionUnbound
#'
#' @description
#' Compute the fraction unbound in vitro
#'
#' @param partitionQSPR type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers & Rowland, Schmidtt)
#' @param inVitroCompartment a list of values describing the in vitro compartment created by `getInVitroCompartment`
#' @param logLipo LogP or LogMA of the compound
#' @param hlcAt Henry's Law Constant in atm/(m3*mol)
#' @param BP Blood plasma ratio
#' @param fu In Vivo Fraction Unbound in plasma from literature
#' @param ionization Vector of length 3 with ionization class, acid, neutral and base, if not input then it is c(0,0,0)
#' @param pKa vector of pkA of the compound
#' @param cCells concentration of cells million /mL
#' @param cMicro concentration of microsomes mg/mL
#' @param volMedium
#'
#' @return  fuInvitro and possible warning for evaporation
#' @export
#'
#' @examples
#' FractionUnbound(
#'   partitionQSPR = "All PK-Sim Standard", logLipo = 3, ionization = c("acid", 0, 0),
#'   typeSystem = "hepatocytes", FBS = 0, microplateType = 96,
#'   volMedium = 0.22, pKa = c(6, 0, 0), hlcAt = 1E-6, cCells = 2
#' )
#'
#' @details
#'
#' mayeb consider to have average data..
FractionUnbound <- function(partitionQSPR, logLipo, ionization,
                            typeSystem, FBS, microplateType,
                            volMedium, pKa = NULL, hlcAt = NULL, fu = NULL, BP = NULL,
                            cMicro = NULL, cCells = NULL) {
  # Have default value of medium volume considering the proportion of the common sim.
  # Have default value of fu
  # How for BP

  # check if the arguments are valid
  rlang::arg_match(partitionQSPR, c(
    "All Poulin and Theil",
    "Poulin and Theil + fu",
    "All Berezhkovskiy",
    "Berezhkovskiy + fu",
    "All PK-Sim Standard",
    "PK-Sim Standard + fu",
    "Rodgers & Rowland + fu",
    "All Schmitt",
    "Schmitt + fu"
  ))

  # run function to get ionization factors
  fneutral <- getIonization(ionization, pKa)
  X <- fneutral["fneutral_plasma"] # Interstitial tissue
  Y <- fneutral["fneutral_cells"] # intracellular

  kPro <- (0.81 + 0.11 * 10^logLipo) / 24.92 * 5.0

  # Calculate air-water partition coefficient
  # Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  # last factor is to convert form atm to Pa
  kAir <- hlcAt / (0.08206 * 310) * 101325

  # Calculate plastic partitioning
  kPlasticFischer <- 10**(logLipo * 0.47 - 4.64)
  kPlasticKramer <- 10**(logLipo * 0.97 - 6.94)
  kPlastic <- as.double(mean(kPlasticFischer, kPlasticKramer) * 1 / (1 + Y))

  # get in vitro compartments----------------------------------------------------
  if (typeSystem == "microsomes") {
    InVitroCompartment <- getInVitroCompartment("microsomes",
      FBS = FBS,
      microplateType = microplateType,
      volMedium = volMedium, cMicro = cMicro
    )
  } else if (typeSystem == "hepatocytes") {
    InVitroCompartment <- getInVitroCompartment("hepatocytes",
      FBS = FBS,
      microplateType = microplateType,
      volMedium = volMedium, cCells = cCells
    )
  }

  cCellNL <- as.double(InVitroCompartment["cCellNL"])
  cCellNPL <- as.double(InVitroCompartment["cCellNPL"])
  cCellAPL <- as.double(InVitroCompartment["cCellAPL"])
  cCellPro <- as.double(InVitroCompartment["cCellPro"])
  cMediumNL <- as.double(InVitroCompartment["cMediumNL"])
  cMediumNPL <- as.double(InVitroCompartment["cMediumNPL"])
  cMediumPro <- as.double(InVitroCompartment["cMediumPro"])
  saPlasticVolMedium <- as.double(InVitroCompartment["saPlasticVolMedium"])
  volAir_L <- as.double(InVitroCompartment["volAir_L"])

  # QSPRs for calculating partitioning in in vitro------------------------------

  if (partitionQSPR == "All Poulin and Theil" || partitionQSPR == "All Berezhkovskiy") {
    # Calculate protein partitioning

    # Calculate lipid partitioning
    kNL <- 10^logLipo
    kPL <- 0.3 * 10^logLipo + 0.7

    fuInvitro <- as.double(1 / (1 + kNL * (cCellNL + cMediumNL) +
      kPL * (cCellNPL + cMediumNPL) +
      kPlastic * saPlasticVolMedium))
  } else if (partitionQSPR == "Poulin and Theil + fu" || partitionQSPR == "Berezhkovskiy + fu") {
    # Calculate lipid partitioning
    kNL <- 10^logLipo
    kNPL <- 0.3 * 10^logLipo + 0.7

    fuInvitro <- 1 / (1 + kNL * cCellNL
      + kNPL * cCellNPL
      + kPlastic * saPlasticVolMedium
      + (1 / fu - 1) * FBS)
  } else if (partitionQSPR == "All PK-Sim Standard") {
    kNL <- 10^logLipo

    # assume all neutral lipids have same binding
    fuInvitro <- as.double(1 / (1 + kNL * (cCellNL + cMediumNL + cCellNPL + cMediumNPL)
      + kPlastic * saPlasticVolMedium
      + kPro * (cCellPro + cMediumPro)))
  } else if (partitionQSPR == "PK-Sim Standard + fu") {
    kNL <- 10^logLipo


    fuInvitro <- as.double(1 / (1 + kNL * (cCellNL + cCellNPL)
      + kPlastic * saPlasticVolMedium
      + kPro * cCellPro
      + (1 / fu - 1) * FBS))
  } else if (partitionQSPR == "Rodgers & Rowland + fu") {
    # RR can only be used by using fu
    # partition into acid phospholipids is only considered if chemical is a strong base
    # this is done by considering X as 0 for acids and neutral chemic.

    kOW <- 10^logLipo
    kNL <- kOW * (1 / (1 + Y))
    Hema <- 0.45
    kpuBC <- (Hema - 1 + BP) / (Hema * fu)
    fiwBC <- 0.63
    fnlBC <- 0.003
    fnpBC <- 0.0059
    APbc <- 0.57 # acidic phspholipids in blood cells

    KAPL_1 <- max(0, kpuBC -
      (1 + Y) / (1 + X) * fiwBC -
      (kNL * fnlBC + (0.3 * kNL + 0.7) * fnpBC))

    kAPL <- KAPL_1 * (1 + Y) / APbc / Y

    fuInvitro <- as.double(1 / (1 + kNL * (cCellNL + cMediumNL) +
      (kNL * 0.3 + 0.7) * (cCellNPL + cMediumNPL) +
      kAPL * (cCellAPL) * X / (1 + X) +
      kPlastic * saPlasticVolMedium))
  } else if (partitionQSPR == "All Schmitt") {
    LogP <- logLipo
    ionParamSchmitt <- getIonizationSchmitt(ionization, pKa)
    logD_Factor <- ionParamSchmitt["logD_Factor"]
    kAPLpHFactor <- ionParamSchmitt["kAPLpHFactor"]
    LogD <- LogP + log10(logD_Factor)
    kNL <- 10**LogD
    kNPL <- 10**LogP
    kAPL <- kNPL * kAPLpHFactor

    fuInvitro <- as.double(1 / (1 + kNL * (cCellNL + cMediumNL) +
      kNPL * (cCellNPL + cMediumNPL) +
      kAPL * (cCellAPL) +
      kPro * (cCellPro + cMediumPro)))
  } else if (partitionQSPR == "Schmitt + fu") {
    LogP <- logLipo

    ionParamSchmitt <- getIonizationSchmitt(ionParam, pKa)
    logD_Factor <- ionParamSchmitt["logD_Factor"]
    kAPLpHFactor <- ionParamSchmitt["kAPLpHFactor"]
    LogD <- LogP + log10(logD_Factor)
    kNPL <- 10**LogP
    kNL <- 10**LogD
    kAPL <- kNPL * kAPLpHFactor

    fuInvitro <- as.double(1 / (1 + kNL * cCellNL +
      kNPL * cCellNPL +
      kAPL * cCellAPL +
      kPro * cCellPro +
      (1 / fu - 1) * FBS))
  } else {}


  # Warning for volatility
  volatility(fuInvitro, kAir, volAir_L)

  return(fuInvitro)
}


getIonizationSchmitt <- function(ionization, pKa) {
  pH <- 7.4

  ionParam <- c(0, 0, 0)
  for (i in seq(1, 3)) {
    if (ionization[i] == "acid") {
      ionParam[i] <- 1
    } else if (ionization[i] == "base") {
      ionParam[i] <- -1
    } else {
      ionParam[i] <- 0
    }
  }

  # Calculate the fraction neutral
  # conditional if molecule is neutral
  if (abs(ionParam[1]) == 1) {
    F1 <- 1 / (1 + 10^(ionParam[1] * (pKa[1] - pH)))
  } else {
    F1 <- 1
  }

  if (abs(ionParam[2]) == 1) {
    F2 <- 1 / (1 + 10^(ionParam[2] * (pKa[2] - pH)))
  } else {
    F2 <- 1
  }

  if (abs(ionParam[2]) == 1) {
    F3 <- 1 / (1 + 10^(ionParam[3] * (pKa[3] - pH)))
  } else {
    F3 <- 1
  }

  # fraction neutral
  K1 <- F1 * F2 * F3
  # fraction with one ionized group
  K2 <- (1 - F1) * F2 * F3
  # fraction with one ionized group
  K3 <- F1 * (1 - F2) * F3
  # fraction with one ionized group
  K4 <- F1 * F2 * (1 - F3)
  # fraction ionized with two groups
  K5 <- (1 - F1) * (1 - F2) * F3
  # fraction ionized with two groups
  K6 <- (1 - F1) * F2 * (1 - F3)
  # fraction ionized with two groups
  K7 <- F1 * (1 - F2) * (1 - F3)
  # Fraction fully ionized
  K8 <- (1 - F1) * (1 - F2) * (1 - F3)

  # taken from schmitt paper
  alpha <- 0.001 # ratio of lipophilciity between the neutral and the charged species of a molecule
  # check eq 9 from Schmitt paper
  logD_Factor <- K1 +
    (K2 + K3 + K4) * alpha^1 +
    K5 * alpha^max(ionParam[1] + ionParam[2], -ionParam[1] - ionParam[2]) +
    K6 * alpha^max(ionParam[1] + ionParam[3], -ionParam[1] - ionParam[3]) +
    K7 * alpha^max(ionParam[3] + ionParam[2], -ionParam[3] - ionParam[2]) +
    K8 * alpha^max(ionParam[1] + ionParam[2] + ionParam[3], -ionParam[1] - ionParam[2] - ionParam[3])

  # check equation 17 and 18 of Schmitt paper
  proportFactorAPL <- 20
  kAPLpHFactor <- K1 +
    K2 * proportFactorAPL^ionParam[1] +
    K3 * proportFactorAPL^ionParam[2] +
    K4 * proportFactorAPL^ionParam[3] +
    K5 * proportFactorAPL^(ionParam[1] + ionParam[2]) +
    K6 * proportFactorAPL^(ionParam[1] + ionParam[3]) +
    K7 * proportFactorAPL^(ionParam[3] + ionParam[2]) +
    K8 * proportFactorAPL^(ionParam[1] + ionParam[2] + ionParam[3])

  return(c("logD_Factor" = logD_Factor, "kAPLpHFactor" = kAPLpHFactor))
}


volatility <- function(fuInvitro, kAir, volAir_L) {
  fuAir <- fuInvitro * kAir * volAir_L
  # if more than 5 % of the chemicals is predicted to evaporate
  # a warning is given
  # this is conservative because HLC is usually for 25 C and not 37 C
  # and the system is not closed but semi-open,
  # hence a prediction of 5 % with this model actually underpredicts how evaporation will occur

  if (fuAir > 0.05) {
    warning("Probable evaporation of test compound")
  } else {}
}
