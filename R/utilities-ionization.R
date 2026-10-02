# Ionization factors in plasma (pH 7.4), in cells (pH 7.0) and in blood cells
# (pH 7.22, used for the Rodgers and Rowland blood cell calibration), used to
# calculate the fraction neutral or ionized. Supports a monoprotic acid, a
# monoprotic base and a zwitterion (one acidic and one basic group); any other
# combination is treated as neutral.
# If there are multiple pKas for acidity or basicity, use the lower value.
# pKb is not the same as pKa: pKa = 14 - pKb.
.calculate_ionization_factors <- function(ionization, pka) {
  groups <- .order_ionizable_groups(ionization, pka)
  ionization <- groups$ionization
  pka <- groups$pka
  # confirm##################
  pH <- 7.4
  pH_cell <- 7.0 # intracellular, average from literature
  pH_blood_cell <- 7.22 # blood cells

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

  ion_factor_at_pH <- function(pH) {
    if (identical(ionParam, c(-1, 0))) {
      # Monoprotic base
      10^(pka[1] - pH)
    } else if (identical(ionParam, c(-1, 1)) | identical(ionParam, c(1, -1))) {
      # monoproticBaseMonoproticAcid
      10^(pka[which(ionParam %in% -1)] - pH) +
        10^(pH - pka[which(ionParam %in% 1)])
    } else if (identical(ionParam, c(1, 0))) {
      # monoprotic acid
      10^(pH - pka[1])
    } else {
      0
    }
  }

  c(
    "ion_factor_plasma" = ion_factor_at_pH(pH),
    "ion_factor_cells" = ion_factor_at_pH(pH_cell),
    "ion_factor_blood_cells" = ion_factor_at_pH(pH_blood_cell)
  )
}

# Put a single ionizable group first, so that c("neutral", "base") with
# pka = c(0, 9) is treated like c("base", "neutral") with pka = c(9, 0).
.order_ionizable_groups <- function(ionization, pka) {
  ionizable <- ionization[1:2] %in% c("acid", "base")
  if (!ionizable[1] && ionizable[2]) {
    ionization <- ionization[c(2, 1)]
    pka <- pka[c(2, 1)]
  }
  list(ionization = ionization, pka = pka)
}
