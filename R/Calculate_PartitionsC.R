#' getInVitroFractionUnbound
#'
#' @description
#' Compute the fraction unbound in vitro
#'
#' @param assumptions type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers & Rowland, Schmidtt)
#' @param inVitroCompartment a list of values describing the in vitro compartment created by `getInVitroCompartment`
#' @param logLipo LogP or LogMA of the compound
#' @param hlcAt Henry's Law Constant in atm/(m3*mol)
#' @param MW Molecular Weight of the compound
#' @param BP Blood plasma ratio
#' @param fu In Vivo Fraction Unbound in plasma from literature
#' @param ionParam Vector of length 3 with ionization class, 1=acid, 0=neutral and -1=base, if not input then it is c(0,0,0)
#' @param pKa pkA of the compound
#'
#'
#' @return  fuInvitro and possible warning for evaporation
#' @export
#'
#' @details
#' Hepatocytes Assays: \eqn{fu_{in\ vitro}=\frac{1}{1 + C_{cells} \times 10^{0.4 \times logP - 1.38}}}

getInVitroFractionUnbound <- function(assumptions, inVitroCompartment, logLipo, hlcAt, MW, BP, ionParam, pKa=NULL,fu = NULL) {

   # check if the arguments are valid
  rlang::arg_match(ionization, c("acid", "base", "neutral"))
  rlang::arg_match(assumptions, c("All Poulin and Theil", "All PK-Sim® Standard", "All Rodgers & Rowland", "All Schmidtt"))

  fNeutral <- getIonization(ionization, pH, pKa)
  KProt <- (0.81 + 0.11 * 10^logLipo) / 24.92 * 5.0

  # Calculate air-water partition coefficient
  # Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  kAir <- hlcAt / (0.08206 * 310)

  # Calculate plastic partitioning
  kPlasticFischer <- 10**(logLipo * 0.47 - 4.64)
  kPlasticKramer <- 10**(logLipo  * 0.97 - 6.94)
  kPlastic <- mean(kPlasticFischer, kPlasticKramer) * fNeutral # check

  # QSPRs for calculating partitioning in in vitro

  if ("assumption" == "All Poulin and Theil"||"assumption" == "All Berezhkovskiy") {
    # Calculate protein partitioning

    # Calculate lipid partitioning
    kNLip <- 10^logLipo
    kPLip <- 0.3*10^logLipo+0.7

    fuInvitro <- 1 / (1 + kNLip * (cNLipCell +cNLipFBS) +
                          kPLip* (cPLipCell+ cPLipFBS) +
                          kPlastic * saPlasticVolMedium)

  } else if ("assumption" == "Poulin and Theil + fu"||"assumption" == "Berezhkovskiy + fu") {

    # Calculate lipid partitioning
    kNLip<- 1.115*logLipo-1.35-LogP_ow_pKa
    kNPLip<- 0.3*10^logLipo+0.7

    fuInvitro <- 1 / (1 + kNLip * cNLipCell
                        + kNPLip * cNPLipCell
                        + kPlastic * saPlasticVolMedium
                        + (1/fu-1)*FBS)

  } else if ("assumption" == "All PK-Sim® Standard") {

    kNLip <- 10^logLipo

     #What is the f_lipids
    fuInvitro=1/(1 + kNLip * ( cNLipCell+cNLipFBS+cPLipCell+cPLipFBS)
                    + KProt*(cProtCell+CProtFBS))


  } else if ("assumption" == "PK-Sim® Standard + fu") {

    kNLip <- 10^logLipo

    #What is the f_lipids
    fuInvitro=1 / ( 1 + kNLip * ( cNLipCell + cPLipCell )
                  + KProt * cProtCell
                  + (1/fu-1)*FBS)

  } else if ("assumption" == "Rodgers & Rowland + fu") {
    #RR can only be used by using fu
    #partition into acid phospholipids is only considered if chemical is a strong base
    #this is done by considering x as 0 for acids and neutral chemic.

    ionization=getIonization(ionParam,pH,pKa)
    X= ionization["X"]
    Y= ionization["Y"]
    Z= ionization["Z"]

    kOW=10^LogP
    kNLip=kOW*(1/(1+Y))
    Hema <- 0.45
    kpuBC <- (Hema - 1 + BP) / (Hema * fu)
    fiwBC <-0.63
    fnlBC <-0.003
    fnpBC <-0.0059

    KAPL_1=max(0,kpuBC -
                (1+Z)/(1+Y)* fiwBC -
                (kNLip* fnlBC + (0.3*kNLip+0.7)* fnpBC))

    KAPL=KAPL_1* (1+Y)/ APbc / Z

    fuInvitro <- 1 / (1 + kNLip * cNLip +
                          (kNLip * 0.3 + 0.7) * cNPLip +
                          kAPL*cAPL*X/(1+Y)+
                          kPla * saPlasticVolMedium)

  } else if ("assumption" == "All Schmitt") {

    LogP=logLipo
    ionParamSchmitt <-getIonizationSchmitt(ionParam,pH,pKa)
    logD_Factor <-ionParamSchmitt["logD_Factor"]
    kAPLpHFactor <-ionParamSchmitt["kAPLpHFactor"]
    LogD <- LogP+Log10(logD_Factor)
    KNPL<- 10**LogP
    KNL <- 10**LogD
    kAPL <- kNPL*kAPLpHFactor

    fuInVitro <- 1 / (1 + kNL * ( cNLipCell+cNLipFBS ) +
                          kNPL * ( cNPLipCell+cNPLipFBS) +
                          kAPL * ( cAPLipCell+cAPLipFBS) +
                          kProt * ( cProtCell +cProtFBS ))

  } else if ("assumption" == "Schmitt and fu") {

    LogP=logLipo

    ionParamSchmitt <-getIonizationSchmitt(ionParam,pH,pKa)
    logD_Factor <-ionParamSchmitt["logD_Factor"]
    kAPLpHFactor <-ionParamSchmitt["kAPLpHFactor"]
    LogD <- LogP+Log10(logD_Factor)
    KNPL <- 10**LogP
    KNL <- 10**LogD
    kAPL <- kNPL*kAPLpHFactor

    fuInVitro <- 1 / (1 + KnL * cNLipCell +
                        KnPL * cNPLipCell +
                        KaPL * cAPLipCell +
                        KProt * cProtCell +
                        (1/fuFBS-1)*FBS)
  } else {}


  # Warning for volatility
  volatility(fuInvitro, kAir, volAir_L)


  return(fuInvitro)
}

getIonization <-function (ionization, pH, pKa){

pH_IW=
pH_BC=7.22

  if (identical(ionParam,c(-1,0,0))){

    #Monoprotic base
    X=10^(pKa[1]-pH_IW)
    Y=10^(pKa[1]-pH)
    Z=10^(pKa[1]-pH_BC)

  } else if (identical(ionParam,c(-1,-1,0))){
    #diprotic base

    X=10^(pKa[1]-pH_IW) + 10^(pKa[1]+pKa[2]-2*pH_IW)
    Y=10^(pKa[1]-pH) + 10^(pKa[1]+pKa[2]-2*pH)
    Z=10^(pKa[1]-pH_BC) + 10^(pKa[1]+pKa[2]-2*pH_BC)

  }else if (identical(ionParam,c(-1,-1,-1))){
    #triproticBase

    X=10^(pKa[1]-pH_IW) + 10^(pKa[1]+pKa[2]-2*pH_IW) +
      10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW)

    Y=10^(pKa[1]-pH) + 10^(pKa[1]+pKa[2]-2*pH) +
      10^(pKa[1]+pKa[2]+pKa[3]-3*pH)

    Z=10^(pKa[1]-pH_BC) + 10^(pKa[1]+pKa[2]-2*pH_BC) +
      10^(pKa[1]+pKa[2]+pKa[3]-3*pH_BC)

  }else if (identical(ionParam,c(-1,1,0))|identical(ionParam,c(1,-1,0))){
    #monoproticBaseMonoproticAcid

    X=10^(pKa[which(ionParam %in% -1)]-pH_IW) +
      10^(pH_IW-pKa[which(ionParam %in% 1)])

    Y=10^(pKa[which(ionParam %in% -1)]-pH) +
      10^(pH-pKa[which(ionParam %in% 1)])

    Z=10^(pKa[which(ionParam %in% -1)]-pH_BC) +
      10^(pH_BC-pKa[which(ionParam %in% 1)])

  }else if (identical(ionParam,c(-1,1,1))|identical(ionParam,c(1,-1,1))|identical(ionParam,c(1,1,-1))){
    #monoproticBaseDiproticAcid

    X=10^(pKa[which(ionParam %in% -1)]-pH_IW) +
      10^(pH_IW-min(pKa[which(ionParam %in% 1)])) +
      10^(2*pH_IW-pKa[which(ionParam %in% 1)][1]-pKa[which(ionParam %in% 1)][2])

    Y=10^(pKa[which(ionParam %in% -1)]-pH) +
      10^(pH-min(pKa[which(ionParam %in% 1)])) +
      10^(2*pH-pKa[which(ionParam %in% 1)][1]-pKa[which(ionParam %in% 1)][2])

    Z=10^(pKa[which(ionParam %in% -1)]-pH_BC) +
      10^(pH_BC-min(pKa[which(ionParam %in% 1)])) +
      10^(2*pH_BC-pKa[which(ionParam %in% 1)][1]-pKa[which(ionParam %in% 1)][2])

  }else if (identical(ionParam,c(-1,-1,1))|identical(ionParam,c(1,-1,-1))|identical(ionParam,c(-1,1,-1))){
    #diproticBaseMonoproticAcid

    X=10^(pH_IW-pKa[which(ionParam %in% 1)]) +
      10^(max(pKa[which(ionParam %in% -1)])-pH_IW) +
      10^(pKa[which(ionParam %in% -1)][1]+pKa[which(ionParam %in% -1)][2]-2*pH_IW)

    Y=10^(pH-pKa[which(ionParam %in% 1)]) +
      10^(max(pKa[which(ionParam %in% -1)])-pH) +
      10^(pKa[which(ionParam %in% -1)][1]+pKa[which(ionParam %in% -1)][2]-2*pH)

    Z=10^(pH_BC-pKa[which(ionParam %in% 1)]) +
      10^(max(pKa[which(ionParam %in% -1)])-pH_BC) +
      10^(pKa[which(ionParam %in% -1)][1]+pKa[which(ionParam %in% -1)][2]-2*pH_BC)

  }else if (identical(ionParam,c(1,1,0))){
    #diproticacid

    X=0
    Y=10^(pH-pKa[1]) + 10^(2*pH-pKa[1]-pKa[2])
    Z=1

  }else if (identical(ionParam,c(1,0,0))){
    #monoprotic acid

    X=0
    Y=10^(pH-pKa[1])
    Z=1

  } else {

    X=0
    Y=0
    Z=1
  }
  return (c("X"=X,"Y"=Y,"Z"=Z))
}


getIonizationSchmitt <- function(ionization, pH, pKa) {
  #Calculate the fraction neutral
  F1 <- 1/(1+10^(ionParam[1]*(pKa[1]-pH)))
  F2 <- 1/(1+10^(ionParam[2]*(pKa[2]-pH)))
  F3 <- 1/(1+10^(ionParam[3]*(pKa[3]-pH)))

  # fraction neutral
  K1<-F1*F2*F3
  #fraction with one ionized group
  K2<-(1-F1)*F2*F3
  #fraction with one ionized group
  K3<-F1*(1-F2)*F3
  #fraction with one ionized group
  K4<-F1*F2*(1-F3)
  #fraction ionized with two groups
  K5<-(1-F1)*(1-F2)*F3
  #fraction ionized with two groups
  K6<-(1-F1)*F2*(1-F3)
  #fraction ionized with two groups
  K7<-F1*(1-F2)*(1-F3)
  #Fraction fully ionized
  K8<-(1-F1)*(1-F2)*(1-F3)

  #taken from schmitt paper
  alpha=-0.003 # ratio of lipophilciity between the neutral and the carged species of a molecule
  #check eq 9 from Schmitt paper
  logD_Factor<-K1+
              (K2 +K3 +K4)* alpha ^ 1 +
                K5* alpha ^ max(ionParam[1]+ionParam[2],-ionParam[1]-ionParam[2]) +
                K6* alpha  ^ max(ionParam[1]+ionParam[3],-ionParam[1]-ionParam[3]) +
                K7* alpha  ^ max(ionParam[3]+ionParam[2],-ionParam[3]-ionParam[2]) +
                K8* alpha  ^ max(ionParam[1]+ionParam[2]+ionParam[3],-ionParam[1]-ionParam[2]-ionParam[3])

  #check equation 17 and 18 of Schmitt paper
  proportFactorAPL=20
  kAPLpHFactor<-K1+
                K2* proportFactorAPL ^ ionParam[1] +
                K3* proportFactorAPL ^ ionParam[2]+
                K4* proportFactorAPL ^ ionParam[3] +
                K5* proportFactorAPL ^ (ionParam[1]+ionParam[2]) +
                K6* proportFactorAPL ^ (ionParam[1]+ionParam[3]) +
                K7* proportFactorAPL ^ (ionParam[3]+ionParam[2]) +
                K8* proportFactorAPL ^ (ionParam[1]+ionParam[2]+ionParam[3])

  return(c("logD_Factor"=logD_Factor,"kAPLpHFactor"=kAPLpHFactor))
}


volatility<-function(fuInvitro,kAir,volAir_L){
            fuAir <- fuInvitro * kAir * volAir_L
            #if more than 5 % of the chemicals is predicted to evaporate
            #a warning is given
            # this is conservative because HLC is usually for 25 C and not 37 C
            # and the system is not closed but semi-open,
            #hence a prediction of 5 % with this model actually underpredicts how evaporation will occur

            if (fuAir > 0.05) {
              warning("Probable evaporation of test compound")
            } else {}
            }

