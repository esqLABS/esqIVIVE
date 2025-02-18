
library(ggplot2)

testFuData<-read.csv("tests/test_fu_microsomes.csv") # Based on collection of Poulins papers

QSARs<-c("All Poulin and Theil",
           "All Berezhkovskiy",
           "All PK-Sim Standard",
           #"Rodgers & Rowland + fu",
           "All Schmitt")
####Rodger and Rolwand is giving problems####
source("R/Calculate_Partitions.R")


for (j in seq(1:length(QSARs)))  {
for (i in seq(1:nrow(testFuData))){

  testFuData[i,10+j]=FractionUnbound(partitionQSPR=QSARs[j],logLipo=as.double(testFuData[i,"LogD"]),
                            hlcAt=0.00001,ionization=c(testFuData[i,"Class"],0,0),
                            typeSystem="microsomes",FBS=0,microplateType=96,
                            volMedium=0.22,pKa=c(as.double(testFuData[i,"pKa"]),0,0),
                            BP=1,fu=0.2,
                            cMicro=as.double(testFuData[i,"Cp"]))


}
}
colnames(testFuData)[seq(11,14)]<-c("All_Poulin_and_Theil",
                                           "All_Berezhkovskiy",
                                           "All_PK_Sim_Standard",
                                           #"Rodgers_Rowland_fu",
                                           "All_Schmidtt")
#plot x as experiemtnal and pred and y
ggplot(testFuData,aes(col=Class))+
  #geom_point(aes(x=Obs,y=Pred.Poulin))+
  #geom_point(aes(x=Obs,y=All_PK_Sim_Standard))
geom_point(aes(x=Obs,y=All_Schmidtt))
