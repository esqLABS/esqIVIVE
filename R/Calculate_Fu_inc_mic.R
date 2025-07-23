#' getInVitroFractionUnbound
#'
#' @description
#' Compute the fraction unbound in vitro
#'
#' @param partition_qspr type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers & Rowland, Schmidtt,
#' then from literature, Poulin, Turner, Austin and Halifax. See thevignette for more details)
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param type_system microsomes or hepatocytes
#' @param fetal_bovine_serum_fraction fraction of serum concentration, values can only go from 0-1
#' @param microplate_type number of wells in the microplate
#' @param volume_medium volume of medium in the well (in mL)
#' @param pka vector of length of 2 with pkA of the compound
#' @param henry_law_constant Henry's Law Constant (in atm/(m3*mol))
#' @param fraction_unbound In Vivo Fraction Unbound in plasma from literature
#' @param blood_plasma_ratio Blood plasma ratio, this parameter is needed for Rodgers and Rowland and Poulin method for basic chemicals
#' @param concentration_microsomes concentration of microsomes (in mg/mL)
#' @param concentration_cells concentration of cells (in million cells/mL)
#'
#' @return  fuInvitro and possible warning for evaporation
#' @export
#'
#' @examples
#'calculate_fu_in_vitro(
#'  partition_qspr = "All PK-Sim Standard", log_lipophilicity = 3, ionization = c("acid", 0),
#'  type_system = "hepatocytes", fetal_bovine_serum_fraction = 0, microplate_type = 96,
#'  volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2)
#'
#'calculate_fu_in_vitro(
#'  partition_qspr = "Poulin and Theil + fu", log_lipophilicity = 3, ionization = c("acid", 0),
#'  type_system = "hepatocytes", fetal_bovine_serum_fraction = 0, microplate_type = 96,fraction_unbound=0.01,blood_plasma_ratio=2,
#'  volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2)
#'
#'calculate_fu_in_vitro(
#'  partition_qspr = "All Schmitt", log_lipophilicity = 0.42, ionization = c("acid", 0),
#'  type_system = "microsomes", fetal_bovine_serum_fraction = 0, microplate_type = 96,fraction_unbound=0.2,blood_plasma_ratio=1,
#'  volume_medium = 0.22, pka = c(6, 0), concentration_microsomes = 2)
#'
#' @details
#'
#' mayeb consider to have average data..
calculate_fu_in_vitro <- function(
  partition_qspr,
  log_lipophilicity,
  ionization,
  type_system,
  fetal_bovine_serum_fraction,
  microplate_type,
  volume_medium,
  pka = NULL,
  henry_law_constant = NULL,
  fraction_unbound = NULL,
  blood_plasma_ratio = NULL,
  concentration_microsomes = NULL,
  concentration_cells = NULL
) {
  # Have default value of medium volume considering the proportion of the common sim.
  # Have default value of fu
  # How for BP

  # check if the arguments are valid
  rlang::arg_match(
    partition_qspr,
    c(
      "Austin",
      "Halifax",
      "Turner",
      "Poulin",
      "All Poulin and Theil",
      "Poulin and Theil + fu",
      "All Berezhkovskiy",
      "Berezhkovskiy + fu",
      "All PK-Sim Standard",
      "PK-Sim Standard + fu",
      "Rodgers & Rowland + fu",
      "All Schmitt",
      "Schmitt + fu",
      "All_literature"
    )
  )

  #make warning that micro and hep should match
  if (
    type_system == "hepatocytes" && !is.null(concentration_cells)
  ) {} else if (
    type_system == "microsomes" && !is.null(concentration_microsomes)
  ) {} else {
    warning("input Type system not matching with cCells or cMicro input")
  }

  # run function to get ionization factors

  ionization_factors <- calculate_ionization_factors(ionization, pka)
  ion_factor_plasma <- ionization_factors["ion_factor_plasma"] # Interstitial tissue
  ion_factor_cells <- ionization_factors["ion_factor_cells"] # intracellular

  protein_partition_1 <- 0.73 * 10^log_lipophilicity - 0.39 # from Endo 2012,dx.doi.org/10.1021/es303379y partition to chicken muscle, R2=0.86
  protein_partition_2 <- 0.163 + 0.0221 * 10^log_lipophilicity #from Schmitt 2008, doi:10.1016/j.tiv.2007.09.010
  kPro <- mean(protein_partition_1, protein_partition_2)

  # Calculate air-water partition coefficient
  #default to a low hlc if it is not given
  if (is.null(henry_law_constant)) {
    henry_law_constant <- 0.000001
  } else {}
  # Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  # last factor is to convert form atm to Pa
  kAir <- henry_law_constant /
    (0.08206 * 310) *
    101325

  # Calculate plastic partitioning
  plastic_partition_fischer <- 10**(log_lipophilicity * 0.47 - 4.64)
  plastic_partition_kramer <- 10**(log_lipophilicity * 0.97 - 6.94)
  kPlastic <- as.double(
    mean(plastic_partition_fischer, plastic_partition_kramer) *
      1 /
      (1 + ion_factor_cells)
  )

  # get in vitro compartments----------------------------------------------------
  if (type_system == "microsomes") {
    in_vitro_compartment <- calculate_in_vitro_compartments(
      "microsomes",
      fetal_bovine_serum_fraction = fetal_bovine_serum_fraction,
      microplateType = microplate_type,
      volMedium_mL = volume_medium,
      cMicro_mgml = concentration_microsomes
    )
  } else if (type_system == "hepatocytes") {
    in_vitro_compartment <- calculate_in_vitro_compartments(
      "hepatocytes",
      fetal_bovine_serum_fraction = fetal_bovine_serum_fraction,
      microplateType = microplate_type,
      volMedium_mL = volume_medium,
      cCells_Mml = concentration_cells
    )
  }

  cCellNL <- as.double(in_vitro_compartment["cCellNL_vvmedium"])
  cCellNPL <- as.double(in_vitro_compartment["cCellNPL_vvmedium"])
  cCellAPL <- as.double(in_vitro_compartment["cCellAPL_vvmedium"])
  cCellPro <- as.double(in_vitro_compartment["cCellPro_vvmedium"])
  cMediumNL <- as.double(in_vitro_compartment["cMediumNL_vvmedium"])
  cMediumNPL <- as.double(in_vitro_compartment["cMediumNPL_vvmedium"])
  cMediumPro <- as.double(in_vitro_compartment["cMediumPro_vvmedium"])
  saPlasticVolMedium <- as.double(in_vitro_compartment[
    "saPlasticVolMedium_m2L"
  ])
  volAir_L <- as.double(in_vitro_compartment["volAir_L"])

  # QSPRs for calculating partitioning in in vitro------------------------------

  if (
    partition_qspr == "All Poulin and Theil" ||
      partition_qspr == "All Berezhkovskiy"
  ) {
    # Calculate protein partitioning

    # Calculate lipid partitioning
    kNL <- 10^log_lipophilicity
    kPL <- 0.3 * 10^log_lipophilicity + 0.7

    fuInvitro <- as.double(
      1 /
        (1 +
          kNL * (cCellNL + cMediumNL) +
          kPL * (cCellNPL + cMediumNPL) +
          kPlastic * saPlasticVolMedium)
    )
  } else if (
    partition_qspr == "Poulin and Theil + fu" ||
      partition_qspr == "Berezhkovskiy + fu"
  ) {
    # Calculate lipid partitioning
    kNL <- 10^log_lipophilicity
    kNPL <- 0.3 * 10^log_lipophilicity + 0.7

    fuInvitro <- 1 /
      (1 +
        kNL * cCellNL +
        kNPL * cCellNPL +
        kPlastic * saPlasticVolMedium +
        (1 / fraction_unbound - 1) * fetal_bovine_serum_fraction)
  } else if (partition_qspr == "All PK-Sim Standard") {
    kNL <- 10^log_lipophilicity

    # assume all neutral lipids have same binding
    fuInvitro <- as.double(
      1 /
        (1 +
          kNL * (cCellNL + cMediumNL + cCellNPL + cMediumNPL) +
          kPlastic * saPlasticVolMedium +
          kPro * (cCellPro + cMediumPro))
    )
  } else if (partition_qspr == "PK-Sim Standard + fu") {
    kNL <- 10^log_lipophilicity

    fuInvitro <- as.double(
      1 /
        (1 +
          kNL * (cCellNL + cCellNPL) +
          kPlastic * saPlasticVolMedium +
          kPro * cCellPro +
          (1 / fraction_unbound - 1) * fetal_bovine_serum_fraction)
    )
  } else if (partition_qspr == "Rodgers & Rowland + fu") {
    # RR can only be used by using fu
    # partition into acid phospholipids is only considered if chemical is a strong base
    # this is done by considering X as 0 for acids and neutral chemic.

    kOW <- 10^log_lipophilicity
    kNL <- kOW * (1 / (1 + blood_plasma_ratio))
    Hema <- 0.45
    kpuBC <- (Hema - 1 + blood_plasma_ratio) / (Hema * fraction_unbound)
    fiwBC <- 0.63
    fnlBC <- 0.003
    fnpBC <- 0.0059
    APbc <- 0.57 # acidic phospholipids in blood cells

    KAPL_1 <- max(
      0,
      kpuBC -
        (1 + blood_plasma_ratio) / (1 + X) * fiwBC -
        (kNL * fnlBC + (0.3 * kNL + 0.7) * fnpBC)
    )

    kAPL <- KAPL_1 * (1 + blood_plasma_ratio) / APbc / blood_plasma_ratio

    fuInvitro <- as.double(
      1 /
        (1 +
          kNL * (cCellNL + cMediumNL) +
          (kNL * 0.3 + 0.7) * (cCellNPL + cMediumNPL) +
          kAPL * (cCellAPL) * X / (1 + X) +
          kPlastic * saPlasticVolMedium)
    )
  } else if (partition_qspr == "All Schmitt") {
    LogP <- log_lipophilicity
    ionization_parameters_schmitt <- calculate_ionization_schmitt(
      ionization,
      pka
    )
    logD_Factor <- ionization_parameters_schmitt["logD_Factor"]
    kAPLpHFactor <- ionization_parameters_schmitt["kAPLpHFactor"]
    LogD <- LogP + log10(logD_Factor)
    kNL <- 10**LogD
    kNPL <- 10**LogP
    kAPL <- kNPL * kAPLpHFactor

    fuInvitro <- as.double(
      1 /
        (1 +
          kNL * (cCellNL + cMediumNL) +
          kNPL * (cCellNPL + cMediumNPL) +
          kAPL * (cCellAPL) +
          kPro * (cCellPro + cMediumPro))
    )
  } else if (partition_qspr == "Schmitt + fu") {
    LogP <- log_lipophilicity

    ionization_parameters_schmitt <- calculate_ionization_schmitt(
      ionization,
      pka
    )
    logD_Factor <- ionization_parameters_schmitt["logD_Factor"]
    kAPLpHFactor <- ionization_parameters_schmitt["kAPLpHFactor"]
    LogD <- LogP + log10(logD_Factor)
    kNPL <- 10**LogP
    kNL <- 10**LogD
    kAPL <- kNPL * kAPLpHFactor

    fuInvitro <- as.double(
      1 /
        (1 +
          kNL * cCellNL +
          kNPL * cCellNPL +
          kAPL * cCellAPL +
          kPro * cCellPro +
          (1 / fraction_unbound - 1) * fetal_bovine_serum_fraction)
    )
  } else if (partition_qspr == "Poulin" & type_system == "microsomes") {
    fraction_unbound_in_vitro <- calculate_fu_poulin(
      ionization = ionization,
      pka = pka,
      ion_factor_plasma = ion_factor_cells,
      ion_factor_cells = ion_factor_cells,
      concentration_cell_neutral_lipids = concentration_cell_neutral_lipids,
      log_lipophilicity = log_lipophilicity
    )
  } else if (partition_qspr == "Austin" && type_system == "microsomes") {
    fraction_unbound_in_vitro <- calculate_fu_austin_microsomes(
      ionization = ionization,
      pka = pka,
      ion_factor = ion_factor_plasma,
      log_lipophilicity = log_lipophilicity,
      concentration_microsomes = concentration_microsomes
    )
  } else if (partition_qspr == "Halifax" && type_system == "microsomes") {
    fraction_unbound_in_vitro <- calculate_fu_halifax(
      ionization = ionization,
      pka = pka,
      ion_factor = ion_factor_plasma,
      log_lipophilicity = log_lipophilicity,
      concentration_microsomes = concentration_microsomes
    )
  } else if (partition_qspr == "Turner" && type_system == "microsomes") {
    fraction_unbound_in_vitro <- calculate_fu_turner(
      ionization = ionization,
      pka = pka,
      ion_factor = ion_factor_plasma,
      log_lipophilicity = log_lipophilicity,
      concentration_microsomes = concentration_microsomes
    )
  } else if (
    partition_qspr == "All_literature" && type_system == "microsomes"
  ) {
    fraction_unbound_in_vitro <- mean(c(
      calculate_fu_poulin(
        ionization = ionization,
        pka = pka,
        ion_factor_plasma = ion_factor_cells,
        ion_factor_cells = ion_factor_cells,
        concentration_cell_neutral_lipids = concentration_cell_neutral_lipids,
        log_lipophilicity = log_lipophilicity
      ),
      calculate_fu_austin_microsomes(
        ionization = ionization,
        pka = pka,
        ion_factor = ion_factor_plasma,
        log_lipophilicity = log_lipophilicity,
        concentration_microsomes = concentration_microsomes
      ),
      calculate_fu_halifax(
        ionization = ionization,
        pka = pka,
        ion_factor = ion_factor_plasma,
        log_lipophilicity = log_lipophilicity,
        concentration_microsomes = concentration_microsomes
      ),
      calculate_fu_turner(
        ionization = ionization,
        pka = pka,
        ion_factor = ion_factor_plasma,
        log_lipophilicity = log_lipophilicity,
        concentration_microsomes = concentration_microsomes
      )
    ))
  } else if (partition_qspr == "Austin" && type_system == "hepatocytes") {
    fraction_unbound_in_vitro <- calculate_fu_austin_hepatocytes(
      ionization = ionization,
      pka = pka,
      ion_factor = ion_factor_plasma,
      log_lipophilicity = log_lipophilicity,
      concentration_cells = concentration_cells
    )
  } else if (partition_qspr == "Poulin" && type_system == "hepatocytes") {
    fraction_unbound_in_vitro <- calculate_fu_poulin(
      ionization = ionization,
      pka = pka,
      ion_factor_plasma = ion_factor_cells,
      ion_factor_cells = ion_factor_cells,
      concentration_cell_neutral_lipids = concentration_cell_neutral_lipids,
      log_lipophilicity = log_lipophilicity
    )
  } else if (partition_qspr == "Kilford" && type_system == "hepatocytes") {
    fraction_unbound_in_vitro <- calculate_fu_kilford(
      ionization = ionization,
      pka = pka,
      ion_factor = ion_factor_plasma,
      log_lipophilicity = log_lipophilicity,
      concentration_cells = concentration_cells
    )
  } else if (
    partition_qspr == "All_literature" && type_system == "hepatocytes"
  ) {
    fraction_unbound_in_vitro <- mean(c(
      calculate_fu_kilford(
        ionization = ionization,
        pka = pka,
        ion_factor = ion_factor_plasma,
        log_lipophilicity = log_lipophilicity,
        concentration_cells = concentration_cells
      ),
      calculate_fu_poulin(
        ionization = ionization,
        pka = pka,
        ion_factor_plasma = ion_factor_cells,
        ion_factor_cells = ion_factor_cells,
        concentration_cell_neutral_lipids = concentration_cell_neutral_lipids,
        log_lipophilicity = log_lipophilicity
      ),
      calculate_fu_austin_hepatocytes(
        ionization = ionization,
        pka = pka,
        ion_factor = ion_factor_plasma,
        log_lipophilicity = log_lipophilicity,
        concentration_cells = concentration_cells
      )
    ))
  } else {}

  # Warning for volatility
  calculate_volatility_correction(fuInvitro, kAir, volAir_L)

  return(fuInvitro)
}


calculate_ionization_schmitt <- function(ionization, pka) {
  pH <- 7.4

  ionParam <- c(0, 0, 0)
  for (i in seq(1, 2)) {
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
    F1 <- 1 / (1 + 10^(ionParam[1] * (pka[1] - pH)))
  } else {
    F1 <- 1
  }

  if (abs(ionParam[2]) == 1) {
    F2 <- 1 / (1 + 10^(ionParam[2] * (pka[2] - pH)))
  } else {
    F2 <- 1
  }

  if (abs(ionParam[2]) == 1) {
    F3 <- 1 / (1 + 10^(ionParam[3] * (pka[3] - pH)))
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
    K8 *
      alpha^max(
        ionParam[1] + ionParam[2] + ionParam[3],
        -ionParam[1] - ionParam[2] - ionParam[3]
      )

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


calculate_volatility_correction <- function(fraction_unbound_in_vitro, air_partition_coefficient, volume_air_l) {
  # Check if any parameters are NULL or NA
  if (
    is.null(fraction_unbound_in_vitro) ||
      is.null(air_partition_coefficient) ||
      is.null(volume_air_l) ||
      is.na(fraction_unbound_in_vitro) ||
      is.na(air_partition_coefficient) ||
      is.na(volume_air_l)
  ) {
    return()
  }

  fuAir <- fraction_unbound_in_vitro * air_partition_coefficient * volume_air_l
  # if more than 5 % of the chemicals is predicted to evaporate
  # a warning is given
  # this is conservative because HLC is usually for 25 C and not 37 C
  # and the system is not closed but semi-open,
  # hence a prediction of 5 % with this model actually underpredicts how evaporation will occur

  if (fuAir > 0.05) {
    warning("Probable evaporation of test compound")
  } else {}
}
