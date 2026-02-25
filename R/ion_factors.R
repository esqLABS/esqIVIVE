#' ion_factors
#'
#' @description
#' Function to calculate fraction unionized
#' If there are multiple pKas for acidity just use the lower value
#' If there are multiple pKas for basicity just use the lower value
#' Mind that pKb is not the same as pKa !
#' If you have pKb, just calculate pKa=14-pKb
#'
#' @param ionization vector of length 2 with ionization class, acid, neutral and base
#' @param pka vector of length 2 with pKa values of the compound
#'
#' @return factors that can be used to calculate the fraction neutral or ionized in plasma and intracellularly
#' @export
#'
#' @examples
#' ion_factors(ionization=c("neutral",0),pka<-c(0,0))
#' ion_factors(ionization=c("acid",0),pka<-c(14,0))
#' ion_factors(ionization=c("base","acid"),pka<-c(5,7))

ion_factors <- function(ionization, pka) {
  # confirm##################
  pH <- 7.4
  pH_cell <- 7.22

  # convert type ionization in 1, 0 and -1
  ionParam <- c(0, 0)
  for (i in seq(1, 2)) {
    if (ionization[i] == "acid") {
      ionParam[i] <- 1
    } else if (ionization[i] == "base") {
      ionParam[i] <- -1
    } else {
      ionParam[i] <- 0
    }
  }

  if (identical(ionParam, c(-1, 0))) {
    # Monoprotic base

    X <- 10^(pka[1] - pH)
    Y <- 10^(pka[1] - pH_cell)
  } else if (identical(ionParam, c(-1, 1)) | identical(ionParam, c(1, -1))) {
    # monoproticBaseMonoproticAcid

    X <- 10^(pka[which(ionParam %in% -1)] - pH) +
      10^(pH - pka[which(ionParam %in% 1)])

    Y <- 10^(pka[which(ionParam %in% -1)] - pH_cell) +
      10^(pH_cell - pka[which(ionParam %in% 1)])
  } else if (identical(ionParam, c(1, 0))) {
    # monoprotic acid

    X <- 10^(pH - pka[1])
    Y <- 10^(pH_cell - pka[1])
  } else {
    X <- 0
    Y <- 0
  }
  return(c("ion_factor_plasma" = X, "ion_factor_cells" = Y))
}
