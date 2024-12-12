#Code to transform from affinity to fraction unb
#It describes how accounting for the affinity to albumin ,
#globulin and membrane lipids such as the one sin lipoproteins
#can be used to calculate the Fu in plasma


#tO DO:
#aDD OTHER SPECIES
#Add possibility to have just QSARs based o lipophlicity
#Add possibility to use PP-LFER QSARs..
#make documentation

#running example
stand_Fu<-convertKintoFu("khsa_kgL"=10^4.48,
                         "kglob_kgL"=10^2.16,
                         "kmemlip_LL"=10^3.51)



###- distribution of Fu_plasma predicted by Kalb -###
convertKintoFu<-function(khsa_kgL,kglob_kgL,kmemlip_LL,fhsa_kgL,fglob_kgL,fmemlip_LL){

  #Average fraction in human plasma ( values form Poulin and xx)
  fhsa_kgL=0.041
  fglob_kgL=0.032
  fmemlip_LL=0.0025

  #assuming density of 1.2 g/L for proteins
  fw<-1-fhsa_kgL/1.2-fmemlip_LL-fglob_kgL/1.2

  K_hsa<-khsa_kgL*fhsa_kgL
  K_glob<-kglob_kgL*fglob_kgL
  K_lip<-kmemlip_LL*fmemlip_LL
  Fu_plasma<-as.double(1/(fw+K_hsa+K_lip+K_glob))

  #just to see proportions in each container
  print(c(Fu_plasma,"falb"=K_hsa*Fu_plasma, "fglob"=K_glob*Fu_plasma,
          "flip"=K_lip*Fu_plasma, "fw"=fw*Fu_plasma))

  return (Fu_plasma)
}


