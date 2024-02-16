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
clearance_IVIVE<-function(typeValue,expData,typeSystem,partitionQSPR,logLipo, hlcAt,ionization,
                          FBS,microplateType,volMedium,
                          pKa=NULL,fu = NULL,cMicro=NULL,cCells=NULL,
                          BP=NULL){

  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes",
                                  "microsomes",
                                 "recombinant enzyme"))

  rlang::arg_match(typeValue, c("decay experimentalcurve",
                                "halfLife",
                                "in vitro clearance parameter" ))

  #deal with optional arguments

  if(missing(pKa)) {
    pKa<-c(0,0,0)
  }else { pKa=pKa }

  if(missing(fu)) {
    fu<-1
  }else { fu=fu }

  if(missing(BP)) {
    BP<-1
  }else { BP=BP }


  #Derive the in vitro clearance value
  if (typeValue=="decay experimentalcurve"){

    kcat<-determineClearance(expData)

  } else if (typeValue=="halfLife") {

    halfLife<-expData
    kcat<-0.693/halfLife

  } else if (typeValue=="in vitro clearance parameter") {

    kcat["Kcat",]<-expData

  } else {

    kcat<-0
  }


  #get fu_in vitro
  if (invitro_system=="microsomes"){
  fuInVitro<-as.double(FractionUnbound(partitionQSPR=partitionQSPR,logLipo=logLipo,
                   hlcAt=hlcAt,ionization=ionization,
                   typeSystem="microsomes",FBS=FBS,microplateType=microplateType,
                   volMedium=0.22,pKa=pKa,
                   BP=BP,fu=fu,
                   cMicro=cMicro))

  nLiver<-N_liver=40  # mg/g
  nInvitro<-cMicro    #mg/mL

  } else if (invitro_system=="hepatocytes"){
  fuInVitro<-as.double(FractionUnbound(partitionQSPR=partitionQSPR,logLipo=logLipo,
                                         hlcAt=hlcAt,ionization=ionization,
                                         typeSystem="hepatocytes",FBS=FBS,microplateType=microplateType,
                                         volMedium=0.22,pKa=pKa,
                                         BP=BP,fu=fu,
                                         cCells=cCells))

  nLiver<-139  # million cells/g liver
  nInvitro<-cCells # million cells/mL assay

  }

  #fraction of liver that is cells
  fintcell_liver=0.67
  #Calculate Clearance for unb fraction
  KcatU<-kcat/fu_invitro
  Cl_spe<-KcatU*nLiver*fintcell_liver/nInvitro

  return(Cl_values)
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

  #Plot for evaluating if fit is reasonable
  diag_plot<-ggplot(clear_curve,aes(x=Time_min,y=Concentration_uM))+
    geom_point()+
    labs(title="fit curve")+
    stat_function(fun = function(x) Kcat_function(x,y0=y0,Kcat=coefficients(fitKcat)), colour = "blue")

  print(diag_plot)

  #Make table with fit Kcat
  fit_95conf=confint(fitKcat)
  kcat=as.double(data.frame("Kcat_permin"=coefficients(fitKcat)["Kcat"],fit_95conf[1],fit_95conf[2]))
  names( kcat)=c("Mean","2.5%_CI","95%_CI")

  return( kcat)


}

# test
# expData<-read.csv("tests/Clearance_example.csv",header=TRUE)
