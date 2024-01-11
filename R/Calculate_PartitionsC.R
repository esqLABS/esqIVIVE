#' getInVitroFractionUnbound
#'
#' @description
#' Compute the fraction unbound in vitro
#'
#' @param LogP LogP of the compound
#' @param hlcAt Henry's Law Constant in atm/(m3*mol)
#' @param MW Molecular Weight of the compound
#' @param pKa pkA of the compound
#' @param ionization type of ionization (acid, base, none)
#' @param assumptions type of assumption used (1: Poulin and Theil, 2: PK-Sim® Standard, 3: Rodgers & Rowland, 4: Schmidtt)
#' @param fu fraction unbound
#' @param BC
#'
#' @return
#' @export
#'
#' @details
#' Hepatocytes Assays: \eqn{fu_{in vitro}=\frac{1}{1 + C_{cells} \times 10^{0.4 \times logP - 1.38}}}
getInVitroFractionUnbound <- function(LogP, hlcAt, MW, pKa, ionization, assumptions, fu, BC) {
  # Calculate ionization
  if (ionization == "acid") {
    fNeutral <- 1 / (1 + 10**(pH - pKa))
  } else if (ionization == "base") {
    fNeutral <- 1 / (1 + 10**(pKa - pH))
  } else {
    fNeutral <- 1
  }

  # Partitions to different components
  kOW <- 10**logP

  # Calculate air-water partition coefficient
  # Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  kAir <- hlcAt / (0.08206 * 310)

  # Calculate plastic partitioning
  kPlasticFischer <- 10**(logP * 0.47 - 4.64)
  kPlasticKramer <- 10**(logP * 0.97 - 6.94)
  kPlastic <- mean(kPlasticFischer, kPlasticKramer) * fNeutral

  # QSPRs for calculating partitioning in in vitro
  if ("assumption" == "Poulin and Theil") {
    # Calculate protein partitioning

    # Calculate lipid partitioning
    kNLip <- fNeutral * kOW

    fuInvitro <- 1 / (1 + kNLip * cLip + kPlastic * saPlasticVolMedium)
  } else if ("assumption" == "PK-Sim® Standard") {

  } else if ("assumption" == "Rodgers & Rowland") {
    Hema <- 0.45
    kpuBC <- (Hema - 1 + BP) / (Hema * fu)
    kAPLip_1 <- kpuBC -
      ((1 + 10**(pKa - pHBC)) / (1 + 10**(pKa - pH)) * fiwBC) -
      (kOW * fnlBC + (0.3 * kOW + 0.7) * fnpBC) / (1 + 10**(pKa - pH))

    kAPBC <- kAPLip_1 * (1 + 10**(pKa - pH)) / (0.5 * 10**(pKa - pHBC))
    kAPLip <- (kAPBC * 10**(pKa - pH)) * fNeutral
    kNLip <- fNeutral * kOW
    if (ionization == "base") {
      fuInvitro <- 1 / (1 + kNLip * cNLip +
        (kNLip * 0.3 + 0.7) * cNPLip +
        kAPLip * cAPLip + # need to check units of APL ( rodgers are in mg/g)
        kPla * saPlasticVolMedium)
    } else {
      fuInvitro <- 1 / (1 + kNLip * cNLip +
        (kNLip * 0.3 + 0.7) * cNPLip +
        kPla * saPlasticVolMedium)
    }
  } else if ("assumption" == "Schmitt") {
    K1 <- F1 * F2 * F3
    F1 * F2 * F3 + K2 * base^CT0 + K3 * base^CT1 + K4 * base^CT2 + K5 * base^(CT0 + CT1) + K6 * base^(CT0 + CT2) + K7 * base^(CT2 + CT1) + K8 * base^(CT0 + CT1 + CT2)

    Kprot <- (0.81 + 0.11 * 10^LogMA) / 24.92 * 5.0
    KnPL <- 10**LogMA

    K_a_PL <- K_n_PL * K_a_PL_pH_Factor
    fuInVitro <- 1 / (1 + KnL * f_nl + KnPL * f_np + KaPL * AP_T + KProt * f_proteins)
  } else if ("assumption" == "Berezhkovskiy") {
    fuInvitro <- 1 / (10^LogP * (V_nlp + 0.3 * V_php) + 0.7 * V_php + V_wp / fu)
  } else if ("assumption" == "Berezhkovskiy") {


  } else {}


  # Warning for volatility
  fuAir <- fuInvitro * kAir * volAir_L
  if (fuAir > 0.1) {
    warning("Probable evaporation")
  } else {}

  return(fuInvitro)
}
