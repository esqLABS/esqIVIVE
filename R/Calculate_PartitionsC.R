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
#' @param ionParam Vector of length 3 with ionization class, 1=acid, 0=neutral and -1=base, if not input then it is c(0,0,0)
#' @param pKa pkA of the compound
#'
#'
#' @return  fuInvitro and possible warning for evaporation
#' @export
#'
#' @details
#'
#'
FractionUnbound <- function(partitionQSPR,logLipo, hlcAt,ionization,
                                      typeSystem,FBS,microplateType,
                                      volMedium,
                                      pKa=NULL,fu = NULL,cMicro=NULL,cCells=NULL,
                                      BP=NULL) {

   # check if the arguments are valid
  rlang::arg_match(partitionQSPR, c("All Poulin and Theil",
                                    "Poulin and Theil + fu",
                                    "All Berezhkovskiy",
                                    "Berezhkovskiy + fu",
                                    "All PK-Sim® Standard",
                                    "PK-Sim® Standard + fu",
                                    "Rodgers & Rowland+ fu",
                                    "All Schmidtt",
                                    "Schmitt and fu"))

  #ionization
  pH<- 7.4

  ionParam<- 0
  for (i in seq(1,3)){
    if(ionization[i]=="acid"){
      ionParam[i]<-1
    } else if (ionization[i]=="base"){
      ionParam[i]<--1
    } else {
      ionParam[i]<-0
    }}

  ionization=getIonization(ionization,pH,pKa)
  X= ionization["X"]
  Y= ionization["Y"]
  Z= ionization["Z"]

  KPro <- (0.81 + 0.11 * 10^logLipo) / 24.92 * 5.0

  # Calculate air-water partition coefficient
  # Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  kAir <- hlcAt / (0.08206 * 310)

  # Calculate plastic partitioning
  kPlasticFischer <- 10**(logLipo * 0.47 - 4.64)
  kPlasticKramer <- 10**(logLipo  * 0.97 - 6.94)
  kPlastic <- as.double(mean(kPlasticFischer, kPlasticKramer) * 1/(1+Y))

  #get in vitro compartments----------------------------------------------------
  source("R/in-vitro-compartment.R")
  if (typeSystem=="microsomes"){

  InVitroCompartment<-getInVitroCompartment("microsomes",FBS=FBS,
                                            microplateType=microplateType,
                                            volMedium=volMedium,cMicro=cMicro)

  } else if (typeSystem=="hepatocytes"){

  InVitroCompartment<-getInVitroCompartment("hepatocytes",FBS=FBS,
                                            microplateType=microplateType,
                                            volMedium=volMedium,cCells=cCells)

  }

  cCellNL = as.double(InVitroCompartment["cCellNL"])
  cCellNPL = as.double(InVitroCompartment["cCellNPL"])
  cCellAPL = as.double(InVitroCompartment["cCellAPL"])
  cCellPro = as.double(InVitroCompartment["cCellPro"])
  cMediumNL = as.double(InVitroCompartment["cMediumNL"])
  cMediumNPL = as.double(InVitroCompartment["cMediumNPL"])
  cMediumPro = as.double(InVitroCompartment["cMediumPro"])
  saPlasticVolMedium = as.double(InVitroCompartment["saPlasticVolMedium"])
  volAir_L = as.double(InVitroCompartment["volAir_L"])

  # QSPRs for calculating partitioning in in vitro------------------------------

  if (partitionQSPR == "All Poulin and Theil"||partitionQSPR == "All Berezhkovskiy") {
    # Calculate protein partitioning

    # Calculate lipid partitioning
    kNL <- 10^logLipo
    kPL<- 0.3*10^logLipo+0.7

    fuInvitro <-as.double(1 / (1 + kNL * (cCellNL+cMediumNL) +
                          kPL* (cCellNPL+cMediumNPL) +
                          kPlastic * saPlasticVolMedium))

  } else if (partitionQSPR == "Poulin and Theil + fu"||partitionQSPR == "Berezhkovskiy + fu") {

    # Calculate lipid partitioning
    kNL<- 10^logLipo
    kNPL<- 0.3*10^logLipo+0.7

    fuInvitro <- 1 / (1 + kNL * cCellNL
                        + kNPL * cCellPL
                        + kPlastic * saPlasticVolMedium
                        + (1/fu-1)*FBS)

  } else if (partitionQSPR == "All PK-Sim® Standard") {

    kNL <- 10^logLipo

     #What is the f_lipids
    fuInvitro=as.double(1/(1 + kNL * (cCellNL+cMediumNL+cCellNPL+cMediumNPL)
                    + kPlastic * saPlasticVolMedium
                    + KPro*(cCellPro+cMediumPro)))


  } else if (partitionQSPR == "PK-Sim® Standard + fu") {

    kNL <- 10^logLipo

    #What is the f_lipids
    fuInvitro=as.double(1 / ( 1 + kNL * ( cCellNL + cCellPL)
                  + kPlastic * saPlasticVolMedium
                  + KPro * cCellPro
                  + (1/fu-1)*FBS))

  } else if (partitionQSPR == "Rodgers & Rowland + fu") {
    #RR can only be used by using fu
    #partition into acid phospholipids is only considered if chemical is a strong base
    #this is done by considering x as 0 for acids and neutral chemic.

    ionFactor=getIonization(ionization,pH,pKa)
    X= ionFactor["X"]
    Y= ionFactor["Y"]
    Z= ionFactor["Z"]

    kOW=10^logLipo
    kNL=kOW*(1/(1+Y))
    Hema <- 0.45
    kpuBC <- (Hema - 1 + BP) / (Hema * fu)
    fiwBC <-0.63
    fnlBC <-0.003
    fnpBC <-0.0059
    APbc<-0.57 # acidic phspholipids in blood cells

    KAPL_1=max(0,kpuBC -
                (1+Z)/(1+Y)* fiwBC -
                (kNL* fnlBC + (0.3*kNL+0.7)* fnpBC))

    KAPL=KAPL_1* (1+Y)/ APbc / Z

    fuInvitro <-as.double( 1 / (1 + kNL * (cCellNL+cMediumNL) +
                          (kNL * 0.3 + 0.7) * (cCellNPL+cMediumNPL) +
                          KAPL*(cCellAPL)*X/(1+Y)+
                          kPlastic * saPlasticVolMedium))

  } else if (partitionQSPR == "All Schmitt") {

    LogP=logLipo
    ionParamSchmitt <-getIonizationSchmitt(ionization,pH,pKa)
    logD_Factor <-ionParamSchmitt["logD_Factor"]
    kAPLpHFactor <-ionParamSchmitt["kAPLpHFactor"]
    LogD <- LogP+log10(logD_Factor)
    KNPL<- 10**LogP
    KNL <- 10**LogP
    kAPL <- kNPL*kAPLpHFactor

    fuInVitro <- as.double(1 / (1 + kNL * ( cCellNL+cMediumNL) +
                          kNPL * ( cCellNPL+cMediumNPL) +
                          kAPL * ( cCellAPL) +
                          kPro * ( cCellPro +cMediumPro)))

  } else if (partitionQSPR == "Schmitt and fu") {

    LogP=logLipo

    ionParamSchmitt <-getIonizationSchmitt(ionization,pH,pKa)
    logD_Factor <-ionParamSchmitt["logD_Factor"]
    kAPLpHFactor <-ionParamSchmitt["kAPLpHFactor"]
    LogD <- LogP+log10(logD_Factor)
    kNPL <- 10**LogP
    kNL <- 10**LogD
    kAPL <- kNPL*kAPLpHFactor

    fuInVitro <- as.double(1 / (1 + kNL * cCellNL+
                          KNPL * cCellNPL+
                          KAPL * cCellAPL+
                          KPro * cCellPro +
                           (1/fu-1)*FBS))
  } else {}


  # Warning for volatility
  volatility(fuInvitro, kAir, volAir_L)


  return(fuInvitro)
}



getIonization <-function (ionParam, pH, pKa){
 #confirm##################
pH_IW=7.4
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

  }else if (identical(ionParam,c(-1,1,0))|identical(ionization,c(1,-1,0))){
    #monoproticBaseMonoproticAcid

    X=10^(pKa[which(ionization %in% -1)]-pH_IW) +
      10^(pH_IW-pKa[which(ionization %in% 1)])

    Y=10^(pKa[which(ionization %in% -1)]-pH) +
      10^(pH-pKa[which(ionization %in% 1)])

    Z=10^(pKa[which(ionization %in% -1)]-pH_BC) +
      10^(pH_BC-pKa[which(ionization %in% 1)])

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

