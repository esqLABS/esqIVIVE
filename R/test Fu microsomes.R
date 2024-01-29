testFuData<-read.csv("tests/test_fu_microsomes.csv")

for (i in seq(1:nrow(testFuData))){

  FractionUnbound(partitionQSPR="Poulin and Theil + fu",logLipo=as.double(testFuData[i,"LogP..37C."]),
                            hlcAt=0.00001,ionization=c(testFuData[i,"Class"],0,0),
                            typeSystem="microsomes",FBS=0,microplateType=96,
                            volMedium=0.22,pKa=c(as.double(testFuData[i,"pKa"]),0,0),
                            BP=1,fu=0.2,
                            cMicro=as.double(testFuData[i,"Cp"]))


}
