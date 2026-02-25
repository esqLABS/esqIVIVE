#' @name calculate_fu_turner 
#' @title Turner algorithm for Fu calculation
#' @description The Turner algorithm for calculating Fu in vitro
#' 
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pKa of the compound
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param conc_mic_mgml concentration of microsomes (in mg/mL)
#'
#' @return fuInvitro
#' @examples
#' calculate_fu_turner("acid",3,3,1)
#' 
#' @export

calculate_fu_turner <- function(
  ionization,
  pKa,
  log_lipophilicity,
  conc_mic_mgml
) {
  
  if (ionization[1] == "base" & pKa[1] > 7) {
    fu_invitro <- 1 /
      (1 + conc_mic_mgml * 10^(0.58 * log_lipophilicity - 2.02))
    
  } else if (ionization[1] == "acid" && pKa[1] < 7) {
    fu_invitro <- 1 /
      (1 + conc_mic_mgml * 10^(0.2 * log_lipophilicity - 1.54))
    
  } else {
    fu_invitro <- 1 /
      (1 + conc_mic_mgml * 10^(0.46 * log_lipophilicity - 1.51))
  }
  return(fu_invitro)
}

#' @name calculate_fu_halifax
#' @title Halifax algorithm for Fu calculation
#' @description The Halifax algorithm for calculating Fu in vitro
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pKa of the compound
#' @param logLipo LogP or LogMA of the compound
#' @param cMicro_mgml concentration of microsomes mg/mL
#' @return fuInvitro
#' @export
#' @examples
#' calculate_fu_halifax(c("base",0),c(3,0),3,1)
#' 
calculate_fu_halifax <- function(
  ionization,
  pKa,
  log_lipophilicity,
  conc_mic_mgml
) {
  ion_factor <- ion_factors(ionization, pKa)["ion_factor_plasma"]
  
  if (ionization[1] == "base" & pKa[1] > 7) {
    log_partition <- log10(as.double(1 / (1 + ion_factor) * 10^log_lipophilicity))
  } else {
    log_partition <- log_lipophilicity
  }
  
  fu_invitro <- 1 /
      (1 + conc_mic_mgml *
          10^(0.072 * log_partition^2 + 0.067 * log_partition - 1.126))

  return(fu_invitro)
}

#' @name calculate_fu_austin_microsomes
#' @title Austin algorithm for microsomes Fu calculation
#' @description The Austin algorithm for calculating Fu in vitro for microsomes
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pKa of the compound
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param conc_mic_mgml concentration of microsomes (in mg/mL)
#' @return fuInvitro
#' @export
#' @examples
#' calculate_fu_austin_microsomes(c("base",0),c(8,0),3,1)
#' 
calculate_fu_austin_microsomes <- function(
  ionization,
  pKa,
  log_lipophilicity,
  conc_mic_mgml
) {
  
  ion_factor <- ion_factors(ionization, pKa)["ion_factor_plasma"]
  if (ionization[1] == "base" & pKa[1] > 7) {
    log_partition <- log10(as.double(1 / (1 + ion_factor ) * 10^log_lipophilicity))
  } else {
    log_partition <- log_lipophilicity
  }

  fu_invitro <- 1 /
    (1 + conc_mic_mgml * 10^(0.56 *  log_partition- 1.41))
  return(fu_invitro)
}

#' @name calculate_fu_austin_hepatocytes
#' @title Austin algorithm for hepatocytes Fu calculation
#' @description The Austin algorithm for calculating Fu in vitro for hepatocytes
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pKa of the compound
#' @param ion_factor ionization factor
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param conc_cell_millionml concentration of hepatocytes (in million cells/mL)
#' @return fuInvitro
#' @export
#' @example calculate_fu_austin_hepatocytes(ionization=c("base",0),pKa=c(3,0),log_lipophilicity=3,conc_cell_millionml=0.5)
calculate_fu_austin_hepatocytes <- function(
  ionization,
  pKa,
  log_lipophilicity,
  conc_cell_millionml
) {
  ion_factor <- ion_factors(ionization, pKa)["ion_factor_plasma"]
  if (ionization[1] == "base" & pKa[1] > 7) {
    log_partition <- log10(as.double(1 / (1 + ion_factor ) * 10^log_lipophilicity))
  } else {
    log_partition <- log_lipophilicity
  }

  fu_invitro <- 1 /
    (1 + conc_cell_millionml * 10^(0.4 * log_partition- 1.38))
  return(fu_invitro)
}

#' @name calculate_fu_kilford
#' @title Kilford algorithm for Fu calculation
#' @description The Kilford algorithm for calculating Fu in vitro
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pKa of the compound
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param conc_cell_millionml concentration of hepatocytes (in million cells/mL)
#' @return fuInvitro
#' @export
#' @example calculate_fu_kilford(ionization=c("base",0),pKa=c(3,0),log_lipophilicity=3,conc_cell_millionml=0.5)
calculate_fu_kilford <- function(
  ionization,
  pKa,
  log_lipophilicity,
  conc_cell_millionml
) {
  ion_factor <- ion_factors(ionization, pKa)["ion_factor_plasma"]
  if (ionization[1] == "base" & pKa[1] > 7) {
    log_partition <- log10(as.double(1 / (1 + ion_factor) * 10^log_lipophilicity))
  } else {
    log_partition <- log_lipophilicity
  }
  volume_ratio <- 0.005 * conc_cell_millionml
  fu_invitro <- 1 /
    (1 +
      125 *
        volume_ratio *
        10^(0.072 * log_partition^2 + 0.067 * log_partition - 1.126))

  return(fu_invitro)
}

#' @name calculate_fu_poulin
#' @title Poulin algorithm for Fu calculation
#' @description The Poulin algorithm for calculating Fu in vitro
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pKa of the compound
#' @param ion_factor_plasma ionization factor
#' @param ion_factor_cells intracellular ionization factor
#' @param blood_plasma blood plasma ratio
#' @param fraction_unbound fraction unbound in plasma
#' @param concentration_cell_neutral_lipids neutral lipid concentration
#' @param log_lipophilicity LogP or LogMA of the compound
#' @return fuInvitro
#' @export
#' @examples
#' calculate_fu_poulin(ionization=c("base",0),pKa=c(6,0),blood_plasma=1,fraction_unbound=0.2,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3)
#' calculate_fu_poulin(ionization=c("neutral",0),pKa=c(0,0),concentration_cell_neutral_lipids=0.03,log_lipophilicity=3)
calculate_fu_poulin <- function(
  ionization,
  pKa,
  blood_plasma=NULL,
  fraction_unbound=NULL,
  concentration_cell_neutral_lipids,
  log_lipophilicity
) {
  neutral_lipid_partition <- 10^log_lipophilicity
  ion_factor_plasma <- ion_factors(ionization, pKa)["ion_factor_plasma"]
  ion_factor_cells <- ion_factors(ionization, pKa)["ion_factor_cells"]
  
  if (ionization[1] == "base" & pKa[1] > 7) {
    fraction_neutral_lipids_erythrocytes <- 0.0024
    fraction_acidic_phospholipids_erythrocytes <- 0.00057
    fraction_water_erythrocytes <- 0.63
    partition_erythrocytes_albumin <- (blood_plasma - (1 - 0.45)) /
      0.45 /
      fraction_unbound
    acidic_phospholipid_partition <- (partition_erythrocytes_albumin -
      ((1 + ion_factor_cells) *
        fraction_water_erythrocytes +
        neutral_lipid_partition * fraction_neutral_lipids_erythrocytes) /
        (1 + ion_factor_plasma)) *
      ((1 + ion_factor_plasma) /
        (ion_factor_cells * fraction_acidic_phospholipids_erythrocytes))

    fu_invitro <- as.double(
      1 /(1 +
          ((neutral_lipid_partition *
            concentration_cell_neutral_lipids +
            ion_factor_plasma * acidic_phospholipid_partition * cCellAPL) /
            (1 + ion_factor_plasma)))
    )
  } else {
    fu_invitro <- as.double(
      1 /
        (1 +
          ((neutral_lipid_partition * concentration_cell_neutral_lipids) /
            (1 + ion_factor_plasma)))
    )
  }
  return(fu_invitro)
}
