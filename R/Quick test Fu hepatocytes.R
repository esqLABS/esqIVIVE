source("R/Calculate_Partitions.R")

FractionUnbound(partitionQSPR="Schmitt + fu",logLipo=3,
                hlcAt=0.00001,ionization=c(0,0,0),
                typeSystem="hepatocytes",FBS=0,microplateType=96,
                volMedium=0.22,pKa=c(0,0,0),
                BP=1,fu=0.2,
                cCells=2)
