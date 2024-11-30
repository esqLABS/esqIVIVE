#'
#' @description
#' IVIVE for clearance
#'
#'
#' @param typeValue describe what type of in vitro data for clearance it is
#' @param units this are the units of the value
#' @param typeSystem microsomes or hepatocytes
#' @param expData depending on the type of value this can be a csv file or just one value, units are min or min-1
#' @param partitionQSPR type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers & Rowland, Schmidtt)
#' @param logLipo LogP or LogMA of the compound
#' @param FBS fraction of FBS or FCS used
#' @param microplateType type of microplate, 96, 24, 48-wells, etc.
#' @param ionization Vector of length 3 with ionization class, acid, neutral and base, if not input then it is c(0,0,0)
#' @param volMedium volume of medium in the well (in mL)
#' @param REF relative expression or activity factor
#' @param hlcAt Henry's Law Constant in atm/(m3*mol)
#' @param pKa pkA of the compound as a vector of length 3
#' @param fu Fraction unbound in plasma
#' @param fu_hep Fraction unbound in the in vitro hepatic system. if not known code will calculate it
#' @param cMicro concentration of microsome protein mg/mL
#' @param cCells concentration of hepatocytes used million cells/mL
#' @param BP Blood plasma ratio
#' @param species values can human and rat for now, defaulting to human
#'
#'@return  Specific clearance parameter to plug in PK-Sim
#'@examples


clearance_IVIVE<-function(typeValue,units,expData,typeSystem,partitionQSPR,logLipo,
                          FBS,microplateType=NULL,ionization=NULL,volMedium=NULL,
                          REF=NULL,hlcAt=NULL,
                          pKa=NULL,fu=NULL,fu_hep=NULL,cMicro=NULL,cCells=NULL,
                          BP=NULL,species=NULL){

  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes",
                                  "microsomes"))

  rlang::arg_match(typeValue, c("decay experimentalcurve",
                                "halfLife",
                                "in vitro clearance parameter" ))

  #defaults for empty variable
  if (is.null(pKa)){
    pKa<-c(0,0,0)
  } else {}

  if (is.null(hlcAt))  {
   hlcAt<-1E-6
  }else{}

  if (is.null(fu))  {
    fu<-0.2
  }else{}

  if (is.null(BP))  {
    BP<-1
  }else{}

  if (is.null(REF)) {
    REF<-1
  }else{}

  if (is.null(species))  {
    species<-"human"
  }else{}

  #make it accordngly to microplate
  if (is.null(volMedium))  {
    volMedium<-1
  }else{}


  if (is.null(fu_hep))  {

      source("R/Calculate_Partitions.R")
      #get fu_hep
      if (typeSystem=="microsomes"){
        fu_hep<-as.double(FractionUnbound(partitionQSPR=partitionQSPR,
                                             logLipo=logLipo,
                                             hlcAt=hlcAt,ionization=ionization,
                                             typeSystem="microsomes",FBS=FBS,
                                             microplateType=microplateType,
                                             volMedium=volMedium,pKa=pKa,
                                             BP=BP,fu=fu, cMicro=cMicro))

        } else if (typeSystem=="hepatocytes"){
        fu_hep<-as.double(FractionUnbound(partitionQSPR=partitionQSPR,
                                             logLipo=logLipo,
                                             hlcAt=hlcAt,ionization=ionization,
                                             typeSystem="hepatocytes",FBS=FBS,
                                             microplateType=microplateType,
                                             volMedium=volMedium,pKa=pKa,
                                             BP=BP,fu=fu,
                                             cCells=cCells))

      } else {}
  } else {}

  if (fu_hep==0) {print("problem fu_hep=0")} else{}

  #Get the scaling factors
  if (typeSystem=="microsomes"){

    nLiver<-40  # mg/g liver
    cInvitro<-cMicro    #mg/mL

  } else if (typeSystem=="hepatocytes"){

    nLiver<-139  # million cells/g liver
    cInvitro<-cCells # million cells/mL assay

  } else {}

  #Derive the in vitro clearance value
  if (typeValue=="decay experimentalcurve"){

    kcat<-determineClearance(expData)[1]
    clinvitro<-kcat/cInvitro

  } else if (typeValue=="halfLife") {

     halfLife<-expData

    #make matrix of calculations depedning on the units
    matrixCalcHalf<-cbind(c("minutes","hours","seconds"),
                          c(1,1/60,60))

    multFactorHalf<-as.double(matrixCalcHalf[which(matrixCalcHalf[,1]==units),2])
    kcat<-0.693/halfLife*multFactorHalf
    clinvitro<-kcat/cInvitro

  } else if (typeValue=="in vitro clearance parameter") {

    #make matrix of calculations depending on the units
    matrixCalClear<-cbind(c("/minutes","/hours","/seconds",
                            "mL/minutes/Millioncells", "uL/minutes/Millioncells", "L/minutes/Millioncells",
                            "mL/hours/Millioncells", "uL/hours/Millioncells", "L/hours/Millioncells",
                            "mL/seconds/Millioncells","uL/seconds/Millioncells","L/seconds/Millioncells",
                            "mL/minutes/cell","uL/minutes/cell",
                            "mL/hours/cell", "uL/hours/cell",
                            "mL/seconds/cell","uL/seconds/cell",
                            "mL/minutes/mg protein", "uL/minutes/mg protein","L/minutes/mg protein",
                            "mL/hours/mg protein", "uL/hours/mg protein","L/hours/mg protein",
                            "mL/seconds/mg protein","uL/seconds/mg protein","L/seconds/mg protein",
                            "mL/minutes","uL/minutes",#not sure if correct
                            "mL/seconds","uL/seconds",
                            "mL/hours","uL/hours"),
                          c(1/cInvitro,1/cInvitro*60,1/cInvitro/60,
                            1,1/1000,1000,
                            1/60,1/60000,1000/60,
                            60,60/1000,60000,
                            1E6,1000,
                            1E6/60,1E6/60000,
                            60*1E6,60*1E6/1000,
                            1,1/1000,1000,
                            1/60,1/60000,1000/60,
                            60,60/1000,60000,
                            1/(cInvitro*volMedium),1/(cInvitro*volMedium)/1000,
                            1/(cInvitro*volMedium)*60,1/(cInvitro*volMedium)*60/1000,
                            1/(cInvitro*volMedium)/60,1/(cInvitro*volMedium)/60000
                            ))

    multFactorClear<-as.double(matrixCalClear[which(matrixCalClear[,1]==units),2])
    clinvitro<-expData*multFactorClear  #mL/min/million cells

  } else {
    clinvitro<-0
  }

  #fraction of liver that is cells
  fintcell_liver=0.67
  #Calculate Clearance for unb fraction
  clinvitroU<- clinvitro/fu_hep
  ClspePermin<-as.double(clinvitroU*nLiver/fintcell_liver*REF)
  print(ClspePermin)

}

#function to determine clearance from experimental curve
determineClearance<-function(expData)  {

  #Load the depletion curve
  clear_curve_xy<-expData
  colnames(clear_curve_xy)<-c("x","y")

  #find the starting concentration
  y0=mean(clear_curve_xy$y[clear_curve_xy$x==0])

  #create clearance model
  Kcat_function <- function(y0,x,Kcat) {
    y=y0*exp(-Kcat*x)
    return(y)
  }

  #fit model
  fitKcat <- nls(y ~  Kcat_function(y0,x,Kcat),
                 data = clear_curve_xy,
                 start=list(Kcat=0.01),
                 trace=TRUE)


  colnames(clear_curve_xy)<-colnames(expData)
  #Plot for evaluating if fit is reasonable
  diag_plot<-ggplot(clear_curve_xy,aes(x=Time_min,y=Concentration_uM))+
    geom_point()+
    labs(title="fit curve")+
    stat_function(fun = function(x) Kcat_function(x,y0=y0,Kcat=coefficients(fitKcat)), colour = "blue")

  print(diag_plot)

  #Make table with fit Kcat
  fit_95conf=confint(fitKcat)
  kcat=as.double(data.frame("Kcat_permin"=coefficients(fitKcat)["Kcat"],fit_95conf[1],fit_95conf[2]))
  names(kcat)=c("Mean","2.5%_CI","95%_CI")

  return(kcat)

}



