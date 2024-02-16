testFuData<-read.csv("tests/test_fu_microsomes.csv") # Based on collection of Poulins papers

QSARs<-c("All Poulin and Theil",
           "Poulin and Theil + fu",
           "All Berezhkovskiy",
           "Berezhkovskiy + fu",
           "All PK-Sim® Standard",
           "PK-Sim® Standard + fu",
           "Rodgers & Rowland + fu",
           "All Schmitt",
           "Schmitt + fu")

source("R/Calculate_Partitions.R")

#THERE IS SOME PROBLEM WITH THE last 3 qsar, SOMEHOW IT DOES NOT PICK THEM UP
for (j in seq(1:length(QSARs)))  {
for (i in seq(1:nrow(testFuData))){

  testFuData[i,11+j]=FractionUnbound(partitionQSPR=QSARs[j],logLipo=as.double(testFuData[i,"LogD"]),
                            hlcAt=0.00001,ionization=c(testFuData[i,"Class"],0,0),
                            typeSystem="microsomes",FBS=0,microplateType=96,
                            volMedium=0.22,pKa=c(as.double(testFuData[i,"pKa"]),0,0),
                            BP=1,fu=0.2,
                            cMicro=as.double(testFuData[i,"Cp"]))


}
}
colnames(testFuData)[seq(11,19)]<-c("All_Poulin_and_Theil",
                                           "Poulin_and_Theil_fu",
                                           "All_Berezhkovskiy",
                                           "Berezhkovskiy_fu",
                                           "All_PK-Sim_Standard",
                                           "PK-Sim_Standard_fu",
                                           "Rodgers_Rowland_fu",
                                           "All_Schmidtt",
                                           "Schmitt_and_fu")
#plot x as experiemtnal and pred and y
ggplot(testFuData)+
  geom_point(aes(x=Obs,y=Pred.Poulin))
