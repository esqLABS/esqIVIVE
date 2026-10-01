# Ionization factors in plasma (pH 7.4) and in cells (pH 7.22), used to
# calculate the fraction neutral or ionized. Supports a monoprotic acid, a
# monoprotic base and a zwitterion (one acidic and one basic group); any other
# combination is treated as neutral.
# If there are multiple pKas for acidity or basicity, use the lower value.
# pKb is not the same as pKa: pKa = 14 - pKb.
.calculate_ionization_factors <- function(ionization, pka) {
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

  c("ion_factor_plasma" = X, "ion_factor_cells" = Y)
}
