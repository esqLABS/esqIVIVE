#Code to derive in vitro clearance and perform IVIVE
#Main author:Susana Proenca
invitro_system="hepatocytes"   #options are hepatocytes or microsomes
N_invitro=3*10**6
type_data="experimental_data"           #options are 1) experimental data, 2) half-life 3) clearance
logP=2
MW=120
Ionization="None"
pK=0
FBS=0
Vol_medium=0.0001 #L

##Mayeb rule can be that if CI are more than 10% we advise checking for the impact with an uncertainty analysis

IVIVE_Clspe<-function(invitro_system,N_invitro,type_data,logP,MW,Ionization,pKa,FBS,Vol_medium)  {

  if(type_data=="experimental_data"){
    #Fit clearance to experimental data
    #Load the depletion curve
    clear_curve=read.csv("Clearance_example.csv",header=TRUE)

    clear_curve_xy=clear_curve
    colnames(clear_curve_xy)=c("x","y")

    #find the starting concentration
    y0=mean(clear_curve_xy$y[clear_curve_xy$x==0])

    #create clearance model
    Kcat_function <- function(y0,x,Kcat) {
      y=y0*exp(-Kcat*x)
      return(y)
    }

    #fit model
    fitKcat <- nls(y ~  Kcat_fuction(y0,x,Kcat),
                   data = clear_curve_xy,
                   start=list(Kcat=0.01),
                   trace=TRUE)

    #Plot for evaluating if fit is reasonable
    diag_plot<-ggplot(clear_curve,aes(x=Time_min,y=Concentration_uM))+
      geom_point()+
      stat_function(fun = function(x) Kcat_function(x,y0=y0,Kcat=coefficients(fitKcat)), colour = "blue")

    print(diag_plot)

    #Make table with fit Kcat
    fit_95conf=confint(fitKcat)
    Cl_values=data.frame("Kcat_permin"=coefficients(fitKcat)["Kcat"],fit_95conf[1],fit_95conf[2])
    names(Cl_values)=c("Mean","2.5%_CI","95%_CI")

  }else if (type_data=="half_life") {
    #if value is given as half/life
    halflife=50  # min################################input
    Cl_values=0.693/halflife
  }else if (type_data=="residual_fraction") {

  }else if (type_data=="clearance") {

    Cl_values= Input_Cl #permin######################input

  }else if (type_data=="clearance") {
    print("error type data")}

  #Universal calculations----------------------------

  fu_invitro<-0.1  # to remove
  #Calculate Kcat for unb fraction
  Cl_values["KcatU",]=Cl_values["Kcat",]/fu_invitro

  #Calculate the specific clearance, it will be in units of per min
  if (invitro_system=="microsomes"){
    N_liver= N_liver=40  # mg/g
  } else if (invitro_system=="hepatocytes"){
    N_liver=139*10**6  # cell/g
  } else {print("Error in vitro system")}
  fintcell_liver=0.67
  #Final value
  Cl_values["Cl_spe",]=Cl_values["KcatU",]*N_liver*fintcell_liver/N_invitro
  return(Cl_values)

}
