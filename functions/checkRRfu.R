#Rodger and Rowland check bases

checkRR<-function(pKa,logLipo,fu){

ionization<-c("base",0,0)

ionParam<-c(0,0,0)
for (i in seq(1,3)){
  if(ionization[i]=="acid"){
    ionParam [i]<-1
  } else if (ionization[i]=="base"){
    ionParam[i]<--1
  } else {
    ionParam[i]<-0
  }}

fneutral=getIonization(ionParam,pKa)
X= fneutral["X"] #Interstitial tissue
Y= fneutral["Y"] #intracellular
Z= fneutral["Z"] #blood cells

kOW=10^logLipo
kNL=kOW*(1/(1+Y))
Hema <- 0.45
kpuBC <- (Hema - 1 + BP) / (Hema * fu)
fiwBC <-0.63
fnlBC <-0.003
fnpBC <-0.0059
APbc<-0.57 # acidic phspholipids in blood cells

KAPL_1=max(0,kpuBC -
             (1+Z)/(1+Y)* fiwBC -
             (kNL* fnlBC + (0.3*kNL+0.7)* fnpBC))

kAPL=KAPL_1* (1+Y)/ APbc / Z

fuInvitro <-as.double( 1 / ( 1+kNL * (cCellNL+cMediumNL) +
                               (kNL * 0.3 + 0.7) * (cCellNPL+cMediumNPL) +
                               kAPL*(cCellAPL)*X/(1+Y)+
                               kPlastic * saPlasticVolMedium))

return(c(kNL,kAPL,fuInvitro))
}



pKa_vector<-seq(from=6,to=10.25,by=0.4)
logLipo<-seq(from=-1,to=6,by=0.3)
fu<-seq(from=0.01,to=0.9,by=0.05)

kNLVector<-0
for (i in seq(1:18)){

  kNLVector[i]<- checkRR(pKa_vector=i)


}
