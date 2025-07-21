#Import dataset from Riley
library(devtools)
#devtools:install(ESQhtpbpk)
pak::pak("esqLABS/ESQhtpbpk")
library(ESQhtpbpk) # package for HTPBPK
library(readxl) # package to read xlsx files that contains input data
library(dplyr) # package for easy manipulation of data
library(tidyr) # package for easy manipulation of data
library(stringr) # package for string manipulation
# packages for plotting
library(ggplot2)
library(plotly)


#import chemical properties
data_path<-system.file("data","Obach_1999_Clmicrosomes.csv",package="esqIVIVE")
Exp_dataset<-read.csv("inst/data/Obach_1999_Clmicrosomes.csv") # to change directory

##scenarios 
#use or not Fu mic
#use or not empirical factor
#permeability limited or not

#Make function to add hemicals
createCompound <- function(chembl_id,scenarios) {

  compoundIdx <- which(Exp_dataset$Drug == chemical_ID)
  
  # Set new compound with compound name corresponding to chembl_id
  comp <- Compound$new(ID = chemical_ID)
  Exp_dataset$Cl[is.na( Exp_dataset$Cl)] <- 0
    comp$setPropertyValue(
      "Chlorine count", 
      value = min(10, Exp_dataset$Cl[compoundIdx])
    )
    
  
  # set molecular weight
  comp$setPropertyValue(
    "Molecular weight", 
    value = as.double(Exp_dataset$MW[compoundIdx]), 
    unit = "g/mol"
  )
  

  # Set lipophilicity based on scenario
  comp$setPropertyValue(
    "Lipophilicity", 
    value = Exp_dataset$LogP37C, 
    unit = "Log Units"
  )
  
  # Set Fu based on scenario
  comp$setPropertyValue(
    "Fraction unbound", 
    value = Exp_dataset$Fraction.Unbound.in.Plasma..fu.[compoundIdx], 
    unit = "fraction"
  )

  # Consider pka depending on scenario (restrict to 0-14 range)

    comp$setPropertyValue(
      name = "pKa value 0", 
      value = min(max(c(0,Exp_dataset$pKa[compoundIdx])), 14)
    )
    comp$setPropertyValue(name="Compound type 0", 
                          value=Exp_dataset$Ionization)
  
  # Include default GFR fraction of 1 based on scenario

    comp$addProcessProperty(
      processType = "GFR",
      propertyName = "GFR",
      parName = "GFR fraction",
      dimension = "Dimensionless",
      value = 1)
   
    
  #Clearance, convert from in vitro
    Cl_min<-clearance_IVIVE(typeValue="in vitro clearance parameter",typeSystem="microsomes",units="mL/minutes/mg protein",
                                    expData=Exp_dataset$Cl_invitro_mlminmg[compoundIdx],
                                    fu_invitro=Exp_dataset$fu_mic[compoundIdx],
                                    ,cMicro_mgml=Exp_dataset$Microsomal_mgml[compoundIdx]
                              ,empirical_scalar="No")
    
  # Add total clearance (as hepatic clearance) also add matching fu and lipophilicity
    comp$addProcessProperty(
      processType = "Liver Plasma Clearance",
      propertyName = "Specific clearance",
      parName = "Specific clearance",
      value = Cl_min, 
      unit = "1/min"
    )
  

  return(comp)
}



studies <- createStudies(scenario = scenarios[1, ])

#does it get clearance?
# run Predictions for compound plasma in PVB
results <- runPredictions(
  studies,
  numberOfCores = 5,
  outputSelections = c("Organism|PeripheralVenousBlood|**|Plasma (Peripheral Venous Blood)"),
  simulationResolution = c(0, max(InputStudies$`t_end [h]`) * 60, 0.1),
  saveResults = FALSE,
  saveSimulation = FALSE,
  queueSize = 200
)