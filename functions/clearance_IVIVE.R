#'
#' @description
#' IVIVE for clearance
#'
#'
#' @param typeValue describe what type of in vitro data for clearance it is
#' @param units this are the units of the value, does not matter typeValue is "decay experimentalcurve"
#' @param typeSystem microsomes or hepatocytes
#' @param expData depending on the type of value this can be a csv file or just one value
#' @param volMedium volume of medium in the well (in mL)
#' @param REF relative expression or activity factor
#' @param fu_invitro Fraction unbound in the in vitro hepatic system. if not known code will calculate it
#' @param cMicro concentration of microsome protein mg/mL
#' @param cCells concentration of hepatocytes used million cells/mL
#' @param species values can human and rat for now, defaulting to human
#' @param tissue tissue of interest, since scaling factors are dependent on the tissue, will default to liver
#'
#'@return  Specific clearance parameter to plug in PK-Sim
#'@examples
#'
#'expData<-read.csv("tests/Clearance_example.csv",header=TRUE)
#'clearance_IVIVE(typeValue="decay experimentalcurve",typeSystem="hepatocytes",cCells=0.5,units="/minutes",
#'                expData=expData,fu_invitro=0.5)
#'
#'example hepatocytes
#' clearance_IVIVE(typeValue="in vitro clearance parameter",typeSystem="hepatocytes",units="mL/minutes/Millioncells",
#' expData=18.27,fu_invitro=0.5,cCells=0.5)
#'
#'
#' example microsomes
#' clearance_IVIVE(typeValue="in vitro clearance parameter",typeSystem="microsomes",units="mL/seconds",
#'                 expData=18.27,fu_invitro=0.5,cMicro=0.5,volMedium=0.5)
#'
#'

clearance_IVIVE<-function(typeValue,units,expData,typeSystem,fu_invitro,tissue=NULL,species=NULL,volMedium=NULL,
                          REF=NULL,cMicro=NULL,cCells=NULL){

  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes",
                                  "microsomes"))

  rlang::arg_match(typeValue, c("decay experimentalcurve",
                                "halfLife",
                                "in vitro clearance parameter" ))

  #defaults for empty variable
  if (is.null(fu_invitro))  {
    fu_invitro<-1
  }else if (fu_invitro==0){
    print("problem fu_invitro=0")
    } else {}

  if (is.null(REF)) {
    REF<-1
  }else{}

  if (is.null(species))  {
    species<-"human"
  }else{}

  #make it accordingly to microplate
  if (is.null(volMedium))  {
    volMedium<-1
  }else{}

  if (is.null(tissue))  {
    tissue<-"liver"
  }else{}

   #get table values  for species
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

    kcat_min<-determineClearance(expData)[1]
    clinvitro<-kcat_min/cInvitro

  } else if (typeValue=="halfLife") {

     halfLife<-expData

    #make matrix of calculations depedning on the units
    matrixCalcHalf<-cbind(c("minutes","hours","seconds"),
                          c(1,1/60,60))

    multFactorHalf<-as.double(matrixCalcHalf[which(matrixCalcHalf[,1]==units),2])
    kcat_min<-0.693/halfLife*multFactorHalf
    clinvitro<-kcat_min/cInvitro

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

    if (length(multFactorClear)==0) {
      print("Not identified units") } else {}

    clinvitro<-expData*multFactorClear  #mL/min/million cells

  } else {
    clinvitro<-0
  }

  #Fraction intracellular
  fintcell_db<-c("brain"=0.96,"fat"=0.82,"gonads"=0.88,
              "heart"=0.76,"kidney"=0.57,"gut"=0.88,"liver"=0.67,
              "lung"=0.23)

  fintcell<-as.double(fintcell_db[which(names(fintcell_db)==tissue)])

  #Calculate Clearance for unb fraction
  clinvitroU<- clinvitro/fu_invitro
  ClspePermin<-as.double(clinvitroU*nLiver/fintcell*REF)

  return("ClspePermin"=ClspePermin)

}

#function to determine clearance from experimental curve
determineClearance<-function(expData)  {

  library(ggplot2)

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



