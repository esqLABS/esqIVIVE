# 
#' Calculations for Fu plasma
#'
#' @description
#' @function correct_fu_Pearce is to correct the Fu base don Pearce correction factor. 
#' @function convertKintoFu is to calculate the Fu in plasma based on the affinty to the different components: albumin ,
#' globulin and membrane lipids such as the ones in lipoproteins
#' @function QSARs_plasma calculates the affinity to different components in plasma based on QSARs
#'
#' @param klip_LL LogP or LogMA of the compound
#' @param kalb_Lkg partition to albumin
#' @param kglob_Lkg partition to globulin
#' @param kmemlip_LL partition to membrane lipids
#' @param species species to be considered, now there is data for human, rat, dog, monkey, rabbit and mouse.Mind that this changes the composition in serum but the specific affinties to albumin still need to be used
#'
#'
#' @return  fU_plasma
#'
#' @examplesCode 
#' convertKintoFu("kalb_Lkg"=10^4.48,
#'                         "kglob_Lkg"=10^2.16,
#'                          "kmemlip_LL"=10^3.51,
#'                          "klip_LL"=100,
#'                          "species"="human")
#'correct_fu_Pearce(fu=0.2,lip_LL=4)   
#'  
#'                     
#'@details
#' for the protein partition
#' if unit partition coefficient is L/L then K_Lkg=K_LL/density_kgL
#' if unit partition coefficient is L/mol, K_Lkg=K_Lmol/MW_gmol*1000gkg




correct_fu_Pearce<-function(fu,lip_LL){
  fNL_plasma<-7E-3 #fraction neutral lipids in plasma
  fu_corrected<- 1/((10^lip_LL)*fNL_plasma+1/fu)
  return("Fu_plasma" = fu_corrected)
}

  
  
### - distribution of Fu_plasma predicted by Kalb -###
convertKintoFu <- function(kalb_Lkg, kglob_Lkg, kmemlip_LL, klip_LL, species) {
  # Average fraction in human plasma
  # values of protein from paper: Factors Influencing the Use and Interpretation of Animal Models
  # in the Development of Parenteral Drug Delivery Systems

  # values for membrane lipids come form Absorption and lipoprotein transport of sphingomyelin

  species_types <- c("human", "rat", "dog", "monkey", "rabbit", "mouse")
  falb_kgL <- c(0.041, 0.031, 0.027, 0.049, 0.039, 0.033)
  fglob_kgL <- c(0.033, 0.035, 0.063, 0.038, 0.018, 0.0587)
  # I considered the rest of protein was globulin
  # g/mL to mL/mL with a lipid density of 0.9 g/ml
  fmemlip_LL <- c(0.0025, 0.0012, 0.0027, 0.0025, 0.00123, 0.00122) / 0.9
  # cholesterol and TG, only have values from human, rat and dog, other values are standard
  flip_LL <- c(0.00196, 0.00072, 0.00123, 0.001, 0.001, 0.001)

  nr_species <- which(species_types == species)
  # assuming density of 1.2 g/L for proteins
  fw <- 1 - falb_kgL[nr_species] / 1.2 - fmemlip_LL[nr_species] - fglob_kgL[nr_species] / 1.2 - flip_LL[nr_species]

  K_alb <- kalb_Lkg * falb_kgL[nr_species]
  K_glob <- kglob_Lkg * fglob_kgL[nr_species]
  K_memlip <- kmemlip_LL * fmemlip_LL[nr_species]
  K_lip <- klip_LL * flip_LL[nr_species]
  Fu_plasma <- as.double(1 / (fw + K_alb + K_lip + K_glob + K_memlip))

  # just to see proportions in each container
  print(c("falb" = K_alb * Fu_plasma, "fglob" = K_glob * Fu_plasma,
    "flip" = (K_lip + K_memlip) * Fu_plasma
  ))

  return("Fu_plasma" = Fu_plasma)
}
