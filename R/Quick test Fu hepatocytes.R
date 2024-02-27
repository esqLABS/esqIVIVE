source("R/Calculate_Partitions.R")

FractionUnbound(partitionQSPR="All PK-Sim Standard",logLipo=3,ionization=c("acid",0,0),
                typeSystem="hepatocytes",FBS=0,microplateType=96,
                volMedium=0.22,pKa=c(6,0,0),hlcAt=1E-6,
                cCells=2)
