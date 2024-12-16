#Code to transform from affinity to fraction unb
#It describes how accounting for the affinity to albumin ,
#globulin and membrane lipids such as the ones in lipoproteins
#can be used to calculate the Fu in plasma

#Important links
#references to QSARs or databses of Fu.
#fup calculator (https://drumap.nibiohn.go.jp/fup/).

#for the protein partition
#if unit partition coefficient is L/L then K_Lkg=K_LL/density_kgL
#if unit partition coefficient is L/mol, K_Lkg=K_Lmol/MW_gmol*1000gkg

#tO DO:

#Add possibility to use PP-LFER QSARs.-need to check how ionization is considered
#make documentation

#running example
stand_Fu<-convertKintoFu("khsa_kgL"=10^4.48,
                         "kglob_kgL"=10^2.16,
                         "kmemlip_LL"=10^3.51,
                         "klip_LL"=100,
                         "species"="human")


###- distribution of Fu_plasma predicted by Kalb -###
convertKintoFu<-function(khsa_Lkg,kglob_Lkg,kmemlip_LL,klip_LL,species){

  #Average fraction in human plasma
  #values of protein from paper: Factors Influencing the Use and Interpretation of Animal Models
  #in the Development of Parenteral Drug Delivery Systems

  #values for membrane lipids come form Absorption and lipoprotein transport of sphingomyelin

  species_types<-c("human","rat","dog","monkey","rabbit","mouse")
  fhsa_kgL=c(0.041,0.031,0.027,0.049,0.039,0.033)
  fglob_kgL=c(0.033,0.035,0.063,0.038,0.018,0.0587)
  #I considered the rest of protein was globulin
  #g/mL to mL/mL with a lipid density of 0.9 g/ml
  fmemlip_LL=c(0.0025,0.0012,0.0027,0.0025,0.00123,0.00122)/0.9
  #cholesterol and TG, only have values from human, rat and dog, other values are standard
  flip_LL=c(0.00196,0.00072,0.00123,0.001,0.001,0.001)

  nr_species<-which(species_types==species)
  #assuming density of 1.2 g/L for proteins
  fw<-1-fhsa_kgL[nr_species]/1.2-fmemlip_LL[nr_species]-fglob_kgL[nr_species]/1.2-flip_LL[nr_species]

  K_hsa<-khsa_Lkg*fhsa_kgL[nr_species]
  K_glob<-kglob_Lkg*fglob_kgL[nr_species]
  K_memlip<-kmemlip_LL*fmemlip_LL[nr_species]
  K_lip<-klip_LL*flip_LL[nr_species]
  Fu_plasma<-as.double(1/(fw+K_hsa+K_lip+K_glob+K_memlip))

  #just to see proportions in each container
  print(c(Fu_plasma,"falb"=K_hsa*Fu_plasma, "fglob"=K_glob*Fu_plasma,
          "flip"=(K_lip+K_memlip)*Fu_plasma))

  return ("Fu_plasma"=Fu_plasma)
}

###internal QSARs###------------------------------------------------------------
#for neutral chemicals with logP >4 use the logP QSAR
#for acidic phenols, cabroxylic acids, pyridine and amines you can use the PPLFER.

#Examples-
# QSARs_plasma("logP",1.62,"pKa"=c(0,0),"ionization"=c(0,0))

# QSARs_plasma("PPLFER",3,"pKa"=c(2,0),"ionization"=c("acid",0),
#              "LFER_E"=1.22,"LFER_S"=0.86,"LFER_A"=0.61,"LFER_B"=0.09,
#              "LFER_V"=1.39)


QSARs_plasma<-function(QSAR,logP,pKa,ionization,
                       LFER_E=NULL,LFER_B=NULL,LFER_A=NULL,LFER_S=NULL,LFER_V=NULL){
  #run function to get ionization factors
  source("R/Ionization.R")
  fneutral=getIonization(ionization,pKa)
  X= fneutral["fneutral_plasma"] #Interstitial tissue
  Y= fneutral["fneutral_cells"] #intracellular

  if (QSAR=="logP"){

    logD<-logP*1/(1+X)
    kmemlip_LL<-10^logD
    #for albumin we are not correcting for ionization since acid molecules also bind albumin
    khsa_kgL<-0.163+0.0221*10^logP  #Schmitt equation, check where is it based and which equation was used in the VCBA
    kglob_kgL<-0.163+0.0221*10^logD

  }else if (QSAR=="PPLFER"){
    #Add LFER_a
    LFER_Ei=0.15+LFER_E
    LFER_Vi=-0.0215+LFER_V
    LFER_Bi=2.15-0.204*LFER_S+1.217*LFER_B+0.314*LFER_V
    LFER_Si= 1.224+0.908*LFER_E+0.827*LFER_S+0.453*LFER_V
    LFER_Ai=-0.208-0.058*LFER_S+0.0354*LFER_A+0.076*LFER_V
    LFER_J=1.793+0.267*LFER_E-0.195*LFER_S+0.35*LFER_V

    #check units
    #check possibly appli limit, range chemicals...
    kmemlip_LL=10^(0.29+0.74*LFER_E-0.72*LFER_S-3.63*LFER_B+3.3*LFER_V)
    kbsa_kgL=0.29+0.36*LFER_E-0.26*LFER_S-3.23*LFER_B+2.82*LFER_V
    #equation for ions from https://pubs.acs.org/doi/10.1021/acs.est.5b06176
    kbsa_kgL=0.85+0.63*LFER_Ei-0.63*LFER_Si-0.05*LFER_Ai+2.08*LFER_Bi+2.06*LFER_Vi+3.13*LFER_J
    khsa_kgL=kbsa_kgL
    kmus_kg=-0.24+0.68*LFER_E-0.76*LFER_S-2.29*LFER_B+2.51*LFER_V
    kglob_kgL=kmus_kg
  }

  return(c("kmemlip_LL"=kmemlip_LL,"khsa_kgL"=khsa_kgL,
           "kglob_kgL"=kglob_kgL))
}
