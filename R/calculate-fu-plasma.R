#' Correct the fraction unbound in plasma for binding to neutral lipids
#'
#' @description
#' Applies the Pearce correction to a measured fraction unbound in plasma, to
#' account for binding to the neutral lipids of plasma that is not seen in the
#' measurement.
#'
#' @param fu_plasma Measured fraction unbound in plasma.
#' @param lipophilicity Lipophilicity of the compound (log units), as logP or
#'   log membrane affinity.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @return The corrected fraction unbound in plasma, a single number.
#' @export
#' @examples
#' correct_fu_plasma_pearce(fu_plasma = 0.2, lipophilicity = 4)
correct_fu_plasma_pearce <- function(
  fu_plasma,
  lipophilicity,
  verbose = FALSE
) {
  fNL_plasma <- 7E-3 #fraction neutral lipids in plasma
  fu_corrected <- 1 /
    ((10^lipophilicity) * fNL_plasma + 1 / fu_plasma)

  if (verbose) {
    .print_ivive_result(
      "correct_fu_plasma_pearce",
      inputs = list(fu_plasma = fu_plasma, lipophilicity = lipophilicity),
      result = fu_corrected
    )
  }

  fu_corrected
}

#' Calculate the fraction unbound in plasma from partition coefficients
#'
#' @description
#' Calculates the fraction unbound in plasma from the partition coefficients
#' of the compound to albumin, globulins, membrane lipids (such as those of
#' lipoproteins) and neutral lipids, and the plasma composition of the
#' species. The first three partition coefficients can be predicted with
#' [calculate_plasma_partitions()].
#'
#' @param partition_albumin Partition coefficient to albumin (L/kg).
#' @param partition_globulin Partition coefficient to globulins (L/kg).
#' @param partition_membrane_lipids Partition coefficient to membrane lipids
#'   (L/L).
#' @param partition_neutral_lipids Partition coefficient to neutral lipids
#'   (L/L). The value is used as given, not as a log value.
#' @param species Species whose plasma composition is used: `"human"`,
#'   `"rat"`, `"dog"`, `"monkey"`, `"rabbit"` or `"mouse"`. Only the plasma
#'   composition changes with the species: the partition coefficients must be
#'   those of that species.
#' @param verbose If `TRUE`, print the inputs and the result, with the
#'   fractions bound to albumin, globulins and lipids.
#'
#' @return The fraction unbound in plasma, a single number.
#' @export
#' @examples
#' calculate_fu_plasma(
#'   partition_albumin = 10^4.48,
#'   partition_globulin = 10^2.16,
#'   partition_membrane_lipids = 10^3.51,
#'   partition_neutral_lipids = 100,
#'   species = "human"
#' )
calculate_fu_plasma <- function(
  partition_albumin,
  partition_globulin,
  partition_membrane_lipids,
  partition_neutral_lipids,
  species,
  verbose = FALSE
) {
  # Average fraction in human plasma
  # values of protein from paper: Factors Influencing the Use and Interpretation of Animal Models
  # in the Development of Parenteral Drug Delivery Systems

  # values for membrane lipids come form Absorption and lipoprotein transport of sphingomyelin

  species_types <- c("human", "rat", "dog", "monkey", "rabbit", "mouse")
  species <- rlang::arg_match(species, species_types)
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
  K_lip <- partition_neutral_lipids * flip_LL[nr_species]
  Fu_plasma <- as.double(1 / (fw + K_alb + K_lip + K_glob + K_memlip))

  if (verbose) {
    .print_ivive_result(
      "calculate_fu_plasma",
      inputs = list(
        partition_albumin = partition_albumin,
        partition_globulin = partition_globulin,
        partition_membrane_lipids = partition_membrane_lipids,
        partition_neutral_lipids = partition_neutral_lipids,
        species = species
      ),
      result = c(
        fu_plasma = Fu_plasma,
        fraction_albumin = K_alb * Fu_plasma,
        fraction_globulin = K_glob * Fu_plasma,
        fraction_lipids = (K_lip + K_memlip) * Fu_plasma
      )
    )
  }

  Fu_plasma
}
