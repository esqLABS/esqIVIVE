#' Algorithm4Fumic
#'
#' @description
#' The different algorithms from literature to calculate Fu_mic
#'
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param blood_plasma_ratio Blood plasma ratio, this parameter is needed for Rodgers and Rowland and Poulin method for basic chemicals
#' @param fraction_unbound In Vivo Fraction Unbound in plasma from literature
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pka vector of length of 2 with pkA of the compound
#' @param concentration_microsomes concentration of microsomes (in mg/mL)
#'
#' @return  fuInvitro
#' @export
#'
#' @examples
#' Turner("acid",3,100,3,1)
#' Halifax("base",3,100,3,1)
#' Halifax(c("base",0),c(3,0),100,3,1)
#' Austin("base",3,100,3,1)
#' Poulin("base",pka=3,X=100,Y=120,blood_plasma_ratio=1,fraction_unbound=0.2,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3,concentration_microsomes=1)
#' Poulin("acid",pka=3,X=100,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3,concentration_microsomes=1)
#'
#' @name Turner
#' @title Turner algorithm for Fu calculation
#' @description The Turner algorithm for calculating Fu in vitro
#'
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pka vector of length of 2 with pkA of the compound
#' @param ion_factor ionization factor
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param concentration_microsomes concentration of microsomes (in mg/mL)
#'
#' @return fuInvitro
#'
#' @export
calculate_fu_turner <- function(
  ionization,
  pka,
  ion_factor,
  log_lipophilicity,
  concentration_microsomes
) {
  if (ionization[1] == "base" & pka[1] > 7) {
    fraction_unbound_in_vitro <- 1 /
      (1 +
        concentration_microsomes * 10^(0.58 * log_lipophilicity - 2.02))
  } else if (ionization[1] == "acid" && pka[1] < 7) {
    fraction_unbound_in_vitro <- 1 /
      (1 + concentration_microsomes * 10^(0.2 * log_lipophilicity - 1.54))
  } else {
    fraction_unbound_in_vitro <- 1 /
      (1 +
        concentration_microsomes * 10^(0.46 * log_lipophilicity - 1.51))
  }
  return(fraction_unbound_in_vitro)
}

#' @name calculate_fu_halifax
#' @title Halifax algorithm for Fu calculation
#' @description The Halifax algorithm for calculating Fu in vitro
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pkA of the compound
#' @param X ionization factor
#' @param logLipo LogP or LogMA of the compound
#' @param cMicro_mgml concentration of microsomes mg/mL
#' @return fuInvitro
#' @export
#' @examples
#' calculate_fu_halifax("base",3,100,3,1)
#' calculate_fu_halifax(c("base",0),c(3,0),100,3,1)
calculate_fu_halifax <- function(
  ionization,
  pka,
  ion_factor,
  log_lipophilicity,
  concentration_microsomes
) {
  log_distribution <- 1 / (1 + ion_factor) * 10^log_lipophilicity

  if (ionization[1] == "base" & pka[1] > 7) {
    fraction_unbound_in_vitro <- 1 /
      (1 +
        concentration_microsomes *
          10^(0.072 * log_distribution^2 + 0.067 * log_distribution - 1.126))
  } else {
    fraction_unbound_in_vitro <- 1 /
      (1 +
        concentration_microsomes *
          10^(0.072 * log_lipophilicity^2 + 0.067 * log_lipophilicity - 1.126))
  }
  return(fraction_unbound_in_vitro)
}

#' @name calculate_fu_austin_microsomes
#' @title Austin algorithm for microsomes Fu calculation
#' @description The Austin algorithm for calculating Fu in vitro for microsomes
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pka vector of length of 2 with pkA of the compound
#' @param ion_factor ionization factor
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param concentration_microsomes concentration of microsomes (in mg/mL)
#' @return fuInvitro
#' @export
#' @examples
#' calculate_fu_austin_microsomes("base",3,100,3,1)
calculate_fu_austin_microsomes <- function(
  ionization,
  pka,
  ion_factor,
  log_lipophilicity,
  concentration_microsomes
) {
  if (ionization[1] == "base" & pka[1] > 7) {
    log_partition <- 1 / (1 + ion_factor) * 10^log_lipophilicity
  } else {
    log_partition <- log_lipophilicity
  }

  fraction_unbound_in_vitro <- 1 /
    (1 + concentration_microsomes * 10^(0.56 * log_lipophilicity - 1.41))
  return(fraction_unbound_in_vitro)
}

#' @name calculate_fu_austin_hepatocytes
#' @title Austin algorithm for hepatocytes Fu calculation
#' @description The Austin algorithm for calculating Fu in vitro for hepatocytes
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pka vector of length of 2 with pkA of the compound
#' @param ion_factor ionization factor
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param concentration_cells concentration of hepatocytes (in million cells/mL)
#' @return fuInvitro
#' @export
calculate_fu_austin_hepatocytes <- function(
  ionization,
  pka,
  ion_factor,
  log_lipophilicity,
  concentration_cells
) {
  if (ionization[1] == "base" & pka[1] > 7) {
    log_partition <- 1 / (1 + ion_factor) * 10^log_lipophilicity
  } else {
    log_partition <- log_lipophilicity
  }

  fraction_unbound_in_vitro <- 1 /
    (1 + concentration_cells * 10^(0.4 * log_lipophilicity - 1.38))
  return(fraction_unbound_in_vitro)
}

#' @name calculate_fu_kilford
#' @title Kilford algorithm for Fu calculation
#' @description The Kilford algorithm for calculating Fu in vitro
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pka vector of length of 2 with pkA of the compound
#' @param ion_factor ionization factor
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param concentration_cells concentration of hepatocytes (in million cells/mL)
#' @return fuInvitro
#' @export
calculate_fu_kilford <- function(
  ionization,
  pka,
  ion_factor,
  log_lipophilicity,
  concentration_cells
) {
  if (ionization[1] == "base" & pka[1] > 7) {
    log_partition <- 1 / (1 + ion_factor) * 10^log_lipophilicity
  } else {
    log_partition <- log_lipophilicity
  }
  volume_ratio <- 0.005 * concentration_cells
  fraction_unbound_in_vitro <- 1 /
    (1 +
      125 *
        volume_ratio *
        10^(0.072 * log_lipophilicity^2 + 0.067 * log_lipophilicity - 1.126))

  return(fraction_unbound_in_vitro)
}

#' @name calculate_fu_poulin
#' @title Poulin algorithm for Fu calculation
#' @description The Poulin algorithm for calculating Fu in vitro
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pka vector of length of 2 with pkA of the compound
#' @param ion_factor_plasma ionization factor
#' @param ion_factor_cells intracellular ionization factor
#' @param blood_plasma_ratio blood plasma ratio
#' @param fraction_unbound fraction unbound in plasma
#' @param concentration_cell_neutral_lipids neutral lipid concentration
#' @param log_lipophilicity LogP or LogMA of the compound
#' @return fuInvitro
#' @export
#' @examples
#' calculate_fu_poulin("base",pka=3,ion_factor_plasma=100,ion_factor_cells=120,blood_plasma_ratio=1,fraction_unbound=0.2,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3)
#' calculate_fu_poulin("acid",pka=3,ion_factor_plasma=100,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3)
calculate_fu_poulin <- function(
  ionization,
  pka,
  ion_factor_plasma,
  ion_factor_cells,
  blood_plasma_ratio,
  fraction_unbound,
  concentration_cell_neutral_lipids,
  log_lipophilicity
) {
  neutral_lipid_partition <- 10**log_lipophilicity

  if (ionization[1] == "base" & pka[1] > 7) {
    fraction_neutral_lipids_erythrocytes <- 0.0024
    fraction_acidic_phospholipids_erythrocytes <- 0.00057
    fraction_water_erythrocytes <- 0.63
    partition_erythrocytes_albumin <- (blood_plasma_ratio - (1 - 0.45)) /
      0.45 /
      fraction_unbound
    acidic_phospholipid_partition <- (partition_erythrocytes_albumin -
      ((1 + ion_factor_cells) *
        fraction_water_erythrocytes +
        neutral_lipid_partition * fraction_neutral_lipids_erythrocytes) /
        (1 + ion_factor_plasma)) *
      ((1 + ion_factor_plasma) /
        (ion_factor_cells * fraction_acidic_phospholipids_erythrocytes))

    fraction_unbound_in_vitro <- as.double(
      1 /
        (1 +
          ((neutral_lipid_partition *
            concentration_cell_neutral_lipids +
            ion_factor_plasma * acidic_phospholipid_partition * cCellAPL) /
            (1 + ion_factor_plasma)))
    )
  } else {
    fraction_unbound_in_vitro <- as.double(
      1 /
        (1 +
          ((neutral_lipid_partition * concentration_cell_neutral_lipids) /
            (1 + ion_factor_plasma)))
    )
  }
  return(fraction_unbound_in_vitro)
}
