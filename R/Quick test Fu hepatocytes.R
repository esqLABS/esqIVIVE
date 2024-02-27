source("R/Calculate_Partitions.R")

FractionUnbound(partitionQSPR="All PK-Sim Standard",logLipo=3,
                hlcAt=0.00001,ionization=c("acid",0,0),
                typeSystem="hepatocytes",FBS=0,microplateType=96,
                volMedium=0.22,pKa=c(6,0,0),
                cCells=2)
