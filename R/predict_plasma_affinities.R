#' Predict affinity constant to plasma components based on QSARs
#' 
#' @name predict_plasma_affinities
#'
#' @description
#' Collection of QSAR to obtain affinity to the different components in serum, membrane lipids (memlip)
#' albumin (alb) and globulin (glob)
#'
#' @param QSAR type of QSAR, it can be logP based or PPLFER based (still sorting the ionization)
#' @param logP is the lipophilicity as given by logKow
#' @param pka is a vector of length of 2
#' @param ionization is a vector of length of 2 which should indicate if chemical is neutral, acid basic.
#' There are spots for ionization in case chemical is zwitterion
#' @param LFER_E LFER E parameter
#' @param LFER_B LFER B parameter
#' @param LFER_A LFER A parameter
#' @param LFER_S LFER S parameter
#' @param LFER_V abraham volume
#'
#' @return partition_membrane_lipids (in L/L), partition_albumin (in L/kg) and partition_globulin (in L/kg)
#' @details
#'
#' For neutral chemicals with logP >4 use the logP QSAR
#' For acidic phenols, carboxylic acids, pyridine and amines you can use the PPLFER.
#' fup calculator (https://drumap.nibiohn.go.jp/fup/).
#' To Do:
#' PP-LFER QSARs.-need to check how ionization is considered
#' make documentation
#'
#' @examples
#'
#' predict_plasma_affinities(QSAR="logP", logP=2, pKa=c(3,0), ionization=c("acid",0))
#'
#' predict_plasma_affinities(QSAR="PPLFER", logP=2, pKa=c(3,0), ionization=c("acid",0), LFER_E=1, LFER_B=0, LFER_A=1.5, LFER_S=0.8, LFER_V=2)
#' 
#' @export

predict_plasma_affinities <- function(
  QSAR,
  logP,
  pKa,
  ionization,
  LFER_E = NULL,
  LFER_B = NULL,
  LFER_A = NULL,
  LFER_S = NULL,
  LFER_V = NULL
) {
  fneutral = ion_factors(ionization, pKa)

  X = fneutral["ion_factor_plasma"] #Interstitial tissue

  Y = fneutral["ion_factor_cells"] #intracellular

  if (QSAR == "logP") {
    logD <- logP * 1 / (1 + X)

    if (pKa[1] != 0) {
      kmemlip_LL <- 10^logD
    } else {
      #Yu et al  regression
      kmemlip_LL <- 10^(1.294 + 0.304 * LogP)
    }

    #for albumin we are not correcting for ionization since acid molecules also bind albumin

    kalb_Lkg <- 0.163 + 0.0221 * kmemlip_LL #Schmitt equation

    kglob_Lkg_1 <- 0.163 + 0.0221 * kmemlip_LL #Schmitt equation for general tissue protein

    kglob_Lkg_2 <- 10^(0.37 * logD - 0.29) #based on the eq used in the VCBA

    kglob_Lkg <- mean(kglob_Lkg_1, kglob_Lkg_2)
    
  } else if (QSAR == "PPLFER") {
    #Add LFER_a

    LFER_Ei = 0.15 + LFER_E

    LFER_Vi = -0.0215 + LFER_V

    LFER_Bi = 2.15 - 0.204 * LFER_S + 1.217 * LFER_B + 0.314 * LFER_V

    LFER_Si = 1.224 + 0.908 * LFER_E + 0.827 * LFER_S + 0.453 * LFER_V

    LFER_Ai = -0.208 - 0.058 * LFER_S + 0.0354 * LFER_A + 0.076 * LFER_V

    LFER_J = 1.793 + 0.267 * LFER_E - 0.195 * LFER_S + 0.35 * LFER_V

    #check possibly appli limit, range chemicals...
    kmemlip_LL = 10^(0.29 +
      0.74 * LFER_E -
      0.72 * LFER_S -
      3.63 * LFER_B +
      3.3 * LFER_V)
    kbsa_kgL = 0.29 +
      0.36 * LFER_E -
      0.26 * LFER_S -
      3.23 * LFER_B +
      2.82 * LFER_V
    #equation for ions from https://pubs.acs.org/doi/10.1021/acs.est.5b06176
    kbsa_Lkg = 0.85 +
      0.63 * LFER_Ei -
      0.63 * LFER_Si -
      0.05 * LFER_Ai +
      2.08 * LFER_Bi +
      2.06 * LFER_Vi +
      3.13 * LFER_J
    kalb_Lkg = kbsa_Lkg
    kmus_Lkg = -0.24 +
      0.68 * LFER_E -
      0.76 * LFER_S -
      2.29 * LFER_B +
      2.51 * LFER_V
    kglob_Lkg = kmus_Lkg
  }
  return(c(
    "partition_membrane_lipids" = kmemlip_LL,
    "partition_albumin" = kalb_Lkg,
    "partition_globulin" = kglob_Lkg
  ))
}
