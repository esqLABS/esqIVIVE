expData<-read.csv("tests/Clearance_example.csv",header=TRUE)
library(ggplot2)

source("R/clearance_IVIVE.R")

#example for curve
clearance_IVIVE(typeValue="decay experimentalcurve",expData=expData,
                partitionQSPR="All Schmitt",
                logLipo=2,
                ionization=c("neutral",0,0),
                typeSystem="hepatocytes",
                FBS=0,pKa=c(3,0,0),
                hlcAt=1E-6,
                microplateType=96,
                volMedium=0.2,
                cCells=0.2)

#example for half life
clearance_IVIVE(typeValue="halfLife",expData=3,
                units="hours",
                partitionQSPR="All Poulin and Theil",
                logLipo=2,
                ionization=c("neutral",0,0),
                typeSystem="microsomes",
                FBS=0,pKa=c(3,0,0),
                hlcAt=1E-6,
                microplateType=96,
                volMedium=0.2,
                cMicro=0.2)

#example for clearance values
clearance_IVIVE(typeValue="in vitro clearance parameter",expData=3,
                units="uL/minutes/cell",
                partitionQSPR="All Poulin and Theil",
                logLipo=2,
                ionization=c("neutral",0,0),
                typeSystem="hepatocytes",
                FBS=0,pKa=c(3,0,0),
                hlcAt=0.08,
                microplateType=96,
                volMedium=0.2,
                cCells=0.2)

