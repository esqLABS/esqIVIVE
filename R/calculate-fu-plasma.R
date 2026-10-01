#' @name correct_fu_pls_pearce
#' @title Correct Fu based on Pearce correction factor
#' @description Corrects the Fu based on Pearce correction factor for neutral lipids in plasma
#' @param fraction_unbound Fraction unbound in plasma
#' @param log_lipophilicity LogP or LogMA of the compound
#' @param verbose if TRUE, print the inputs and the resulting Fu_plasma
#' @return Corrected Fu_plasma value
#' @export
#' @examples
#' correct_fu_pls_pearce(fraction_unbound=0.2, log_lipophilicity=4)
#'
  correct_fu_pls_pearce <- function(fraction_unbound, log_lipophilicity, verbose = FALSE) {
  fNL_plasma <- 7E-3 #fraction neutral lipids in plasma
  fu_corrected <- 1 /
    ((10^log_lipophilicity) * fNL_plasma + 1 / fraction_unbound)

  if (verbose) {
    .print_ivive_result(
      "correct_fu_pls_pearce",
      inputs = list(fraction_unbound = fraction_unbound, log_lipophilicity = log_lipophilicity),
      result = fu_corrected
    )
  }
  result<-c("Fu_plasma" = fu_corrected)
  return(result)
}

#' @name calculate_fu_pls_from_Ks
#' @title Calculate Fu in plasma based on affinity to different components
#' @description Calculate the Fu in plasma based on the affinity to the different components: albumin,
#' globulin and membrane lipids such as the ones in lipoproteins
#' @param partition_albumin partition to albumin (in L/kg)
#' @param partition_globulin partition to globulin (in L/kg)
#' @param partition_membrane_lipids partition to membrane lipids (in L/L)
#' @param partition_lipids LogP or LogMA of the compound
#' @param species species to be considered, now there is data for human, rat, dog, monkey, rabbit and mouse. Mind that this changes the composition in serum but the specific affinities to albumin still need to be used
#' @param verbose if TRUE, print the inputs and the resulting Fu_plasma
#' @return Fu_plasma value
#' @export
#' @examples
#' calculate_fu_pls_from_Ks("partition_albumin"=10^4.48,
#'                "partition_globulin"=10^2.16,
#'                "partition_membrane_lipids"=10^3.51,
#'                "partition_lipids"=100,
#'                "species"="human")
#'
  calculate_fu_pls_from_Ks <- function(
  partition_albumin,
  partition_globulin,
  partition_membrane_lipids,
  partition_lipids,
  species,
  verbose = FALSE
) {
  # Average fraction in human plasma
  # values of protein from paper: Factors Influencing the Use and Interpretation of Animal Models
  # in the Development of Parenteral Drug Delivery Systems

  # values for membrane lipids come form Absorption and lipoprotein transport of sphingomyelin

  species_types <- c("human", "rat", "dog", "monkey", "rabbit", "mouse")
  rlang::arg_match(species, species_types)
  falb_kgL <- c(0.041, 0.031, 0.027, 0.049, 0.039, 0.033)
  fglob_kgL <- c(0.033, 0.035, 0.063, 0.038, 0.018, 0.0587)
  # I considered the rest of protein was globulin
  # g/mL to mL/mL with a lipid density of 0.9 g/ml
  fmemlip_LL <- c(0.0025, 0.0012, 0.0027, 0.0025, 0.00123, 0.00122) / 0.9
  # cholesterol and TG, only have values from human, rat and dog, other values are standard
  flip_LL <- c(0.00196, 0.00072, 0.00123, 0.001, 0.001, 0.001)

  nr_species <- which(species_types == species)
  # assuming density of 1.2 g/L for proteins
  fw <- 1 -
    falb_kgL[nr_species] / 1.2 -
    fmemlip_LL[nr_species] -
    fglob_kgL[nr_species] / 1.2 -
    flip_LL[nr_species]

  K_alb <- partition_albumin * falb_kgL[nr_species]
  K_glob <- partition_globulin * fglob_kgL[nr_species]
  K_memlip <- partition_membrane_lipids * fmemlip_LL[nr_species]
  K_lip <- partition_lipids * flip_LL[nr_species]
  Fu_plasma <- as.double(1 / (fw + K_alb + K_lip + K_glob + K_memlip))

  if (verbose) {
    .print_ivive_result(
      "calculate_fu_pls_from_Ks",
      inputs = list(partition_albumin = partition_albumin, partition_globulin = partition_globulin, partition_membrane_lipids = partition_membrane_lipids, partition_lipids = partition_lipids, species = species),
      result = c(
        Fu_plasma = Fu_plasma,
        falb = K_alb * Fu_plasma,
        fglob = K_glob * Fu_plasma,
        flip = (K_lip + K_memlip) * Fu_plasma
      )
    )
  }

  result<-c("Fu_plasma" = Fu_plasma)
  return(result)
}
