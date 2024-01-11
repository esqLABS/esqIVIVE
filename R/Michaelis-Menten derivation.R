#Code to derive in vitro michaelis menten clearance and perform IVIVE
#Main author:Susana Proenca
#Needs to be converted in a function


#libraries
library(ggplot2)



#Calculate in vitro intrinsic clearance values---------------------------------

#Load in vitro data
mm_exp_curve=read.csv("Example_michaelis_menten_curve.csv",header=TRUE)
mm_exp_curve_xy=mm_exp_curve
colnames(mm_exp_curve_xy)=c("x","y")

#fit model
fitmm <- nls(y ~ Vmax * x / (Km + x),
             data = mm_exp_curve_xy,
             start=list(Vmax=max(mm_exp_curve_xy$y),Km=mean(mm_exp_curve_xy$x)),
             trace=TRUE)
fit_95conf=confint(fitmodel)

mm_fuction <- function(x){
  y=coefficients(fitmodel)["Vmax"] * x / (coefficients(fitmodel)["Km"] + x)
  return (y)
}

#Plot for diagnosis
ggplot(mm_exp_curve,aes(x=Concentration_uM),y=Velocity(umol_min_mg_protein))+
  geom_point()+
  stat_function(fun = function(x) mm_fuction(x), colour = "blue")


#Add fit values in dataframe for calculations
fitresults_vmax_km=data.frame("Mean"=c(coefficients(fitmodel)["Km"],
                                       coefficients(fitmodel)["Vmax"]),
                              "2.5_percent"=c(fit_95conf["Km",1],
                                              fit_95conf["Vmax",1]),
                              "95%_percent"=c(fit_95conf["Km",2],
                                              fit_95conf["Vmax",2]))

row.names(fitresults_vmax_km)=c("Km_uM","Vmax_umol_min_mgmicro")

#Describe in vitro experiment---------------------------------------------------

system="microsomes"  #hepatocytes or microsomes
c_invitro=1
c_invitro_units="mg/mL"

#Calculate fraction unbound ---------------------------------------------------

fu_invitro=0.33
#Calculate Vmax for unb fraction with units umol per min * mgmicrossome
fitresults_vmax_km["VmaxU",]=fitresults_vmax_km["Vmax",]/fu_invitro

#Calculate in vivo Cl_u--------------------------------------------------------
sf_cyp=108  #pmol_permgmicro, specific of enzyme

Cl_values["KcatU_permin",]=Cl_values["Kcat_permin",]/fu_invitro

#Calculate the specific clearance, it will be in units of per min
if (invitro_system=="microsomes"){
  N_liver= N_liver=40  # mg/g
} else if (invitro_system=="hepatocytes"){
  N_liver=139*10**6  # cell/g
} else {print("Error in vitro system")}

fu_intracell_liver=0.67

#in specific CYPs---

#Calculate Vmax of unb fraction in min-1
fitresults_vmax_km["VmaxU_scaled",]=fitresults_vmax_km["VmaxU",]*sf_cyp

#Calculate clearance (min-1)  from Vmax and Km
CLspec=fitresults_vmax_km["VmaxU",]/fitresults_vmax_km["Km",]*N_liver


