# Code for sorting ionization

# examples to test
# getIonization(ionization=c("neutral",0),pKa<-c(0,0))
# getIonization(ionization=c("acid",0),pKa<-c(14,0))
# getIonization(ionization=c("base","acid"),pKa<-c(5,7))

# If there are multiple pKas for acidity just use the lower value
# If there are multiple pKas for basicity just use the lower value
# Mind that pKb is not the same as pKa !
# If you have pKb, just calculate pKa=14-pKb

getIonization <- function(ionization, pKa) {
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

    X <- 10^(pKa[1] - pH)
    Y <- 10^(pKa[1] - pH_cell)
  } else if (identical(ionParam, c(-1, 1)) | identical(ionParam, c(1, -1))) {
    # monoproticBaseMonoproticAcid

    X <- 10^(pKa[which(ionParam %in% -1)] - pH) +
      10^(pH - pKa[which(ionParam %in% 1)])

    Y <- 10^(pKa[which(ionParam %in% -1)] - pH_cell) +
      10^(pH_cell - pKa[which(ionParam %in% 1)])
  } else if (identical(ionParam, c(1, 0))) {
    # monoprotic acid

    X <- 10^(pH - pKa[1])
    Y <- 10^(pH_cell - pKa[1])
  } else {
    X <- 0
    Y <- 0
  }
  return(c("fneutral_plasma" = X, "fneutral_cells" = Y))
}
