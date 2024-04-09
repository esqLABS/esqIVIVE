#'
#' @description
#' IVIVE for clearance
#'
#' @param typeValue describe what type of in vitro data for clearance it is
#' @param expData pdepedning on the type of value this can be a csv file or just one value, units are min or min-1
#' @param partitionQSPR type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers & Rowland, Schmidtt)
#' @param inVitroCompartment a list of values describing the in vitro compartment created by `getInVitroCompartment`
#' @param logLipo LogP or LogMA of the compound
#' @param hlcAt Henry's Law Constant in atm/(m3*mol)
#' @param BP Blood plasma ratio
#' @param fu In Vivo Fraction Unbound in plasma from literature
#' @param ionParam Vector of length 3 with ionization class, 1=acid, 0=neutral and -1=base, if not input then it is c(0,0,0)
#' @param pKa pkA of the compound
#' @param cMicro concentration of microsome protein mg/mL
#' @param cCells concentration of hepatocytes used million cells/mL
#'
#'
#'@return  Specific clearance parameter to plug in PK-Sim
#'@examples


clearance_IVIVE<-function(typeValue,expData,typeSystem,partitionQSPR,logLipo,ionization,
                          FBS,microplateType,volMedium,units=NULL,hlcAt=NULL,
                          pKa=NULL,fu=NULL,cMicro=NULL,cCells=NULL,
                          BP=NULL){

  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes",
                                  "microsomes",
                                 "recombinant enzyme"))

  rlang::arg_match(typeValue, c("decay experimentalcurve",
                                "halfLife",
                                "in vitro clearance parameter" ))
  #deal with empty variable
  if (exists("pKa")==FALSE){
    pKa<-c(0,0,0)
  } else {}

  if (exists("hlcAt")==FALSE)  {
   hlcAt<-1E-6
  }else{}

  if (exists("fu")==FALSE)  {
    fu<-0.2
  }else{}

  if (exists("BP")==FALSE)  {
    BP<-1
  }else{}


  source("R/Calculate_Partitions.R")
  #get fu_in vitro
  if (typeSystem=="microsomes"){
    fuInVitro<-as.double(FractionUnbound(partitionQSPR=partitionQSPR,
                                         logLipo=logLipo,
                                         hlcAt=hlcAt,ionization=ionization,
                                         typeSystem="microsomes",FBS=FBS,
                                         microplateType=microplateType,
                                         volMedium=volMedium,pKa=pKa,
                                         BP=BP,fu=fu,
                                         cMicro=cMicro))

    nLiver<-40  # mg/g liver
    nInvitro<-cMicro    #mg/mL

  } else if (typeSystem=="hepatocytes"){
    fuInVitro<-as.double(FractionUnbound(partitionQSPR=partitionQSPR,
                                         logLipo=logLipo,
                                         hlcAt=hlcAt,ionization=ionization,
                                         typeSystem="hepatocytes",FBS=FBS,
                                         microplateType=microplateType,
                                         volMedium=volMedium,pKa=pKa,
                                         BP=BP,fu=fu,
                                         cCells=cCells))

    nLiver<-139  # million cells/g liver
    nInvitro<-cCells # million cells/mL assay

  } else {}

  #Derive the in vitro clearance value
  if (typeValue=="decay experimentalcurve"){

    kcat<-determineClearance(expData)[1]
    clinvitro<-kcat/nInvitro

  } else if (typeValue=="halfLife") {

     halfLife<-expData

    #make matrix of calculations depedning on the units
    matrixCalcHalf<-cbind(c("minutes","hours","seconds"),
                          c(1,1/60,60))

    multFactorHalf<-as.double(matrixCalcHalf[which(matrixCalcHalf[,1]==units),2])
    kcat<-0.693/halfLife*multFactorHalf
    clinvitro<-kcat/nInvitro

  } else if (typeValue=="in vitro clearance parameter") {

    #make matrix of calculations depedning on the units
    matrixCalClear<-cbind(c("/minutes","/hours","/seconds",
                            "mL/minutes/Millioncells", "uL/minutes/Millioncells",
                            "mL/hours/Millioncells", "uL/hours/Millioncells",
                            "mL/seconds/Millioncells","uL/seconds/Millioncells",
                            "mL/minutes/cell","uL/minutes/cell",
                            "mL/hours/cell", "uL/hours/cell",
                            "mL/seconds/cell","uL/seconds/cell"),
                          c(1/nInvitro,1/nInvitro*60,1/nInvitro/60,
                            1,1/1000,
                            1/60,1/60000,
                            60,60/1000,
                            1E6,1000,
                            1E6/60,1E6/60000,
                            60*1E6,60*1E6/1000))

    multFactorClear<-as.double(matrixCalClear[which(matrixCalClear[,1]==units),2])
    clinvitro<-expData*multFactorClear  #mL/min/million cells

  } else {
    clinvitro<-0
  }

  #fraction of liver that is cells
  fintcell_liver=0.67
  #Calculate Clearance for unb fraction
  clinvitroU<- clinvitro/fuInVitro
  ClspePermin<-as.double(clinvitroU*nLiver*fintcell_liver)

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
  names( kcat)=c("Mean","2.5%_CI","95%_CI")

  return(kcat)

}

# test

