# Reference values follow from the Henderson-Hasselbalch relations, with
# plasma pH 7.4, intracellular pH 7.0 and blood cell pH 7.22:
#   acid: ion factor = 10^(pH - pKa),  base: ion factor = 10^(pKa - pH)
# A base + acid compound is the sum of both contributions.

test_that(".calculate_ionization_factors: neutral compound", {
  expect_equal(
    .calculate_ionization_factors(ionization = c("neutral", 0), pka = c(0, 0)),
    c(ion_factor_plasma = 0, ion_factor_cells = 0, ion_factor_blood_cells = 0)
  )
})

test_that(".calculate_ionization_factors: monoprotic acid", {
  expect_equal(
    .calculate_ionization_factors(ionization = c("acid", 0), pka = c(7, 0)),
    c(
      ion_factor_plasma = 10^0.4, # 10^(7.4 - 7)
      ion_factor_cells = 1, # 10^(7.0 - 7)
      ion_factor_blood_cells = 10^0.22 # 10^(7.22 - 7)
    ),
    tolerance = 1e-8
  )
})

test_that(".calculate_ionization_factors: monoprotic acid with a very high pKa", {
  # compared on the log10 scale, since the factors themselves are ~1e-7
  expect_equal(
    log10(.calculate_ionization_factors(
      ionization = c("acid", 0),
      pka = c(14, 0)
    )),
    c(
      ion_factor_plasma = 7.4 - 14,
      ion_factor_cells = 7.0 - 14,
      ion_factor_blood_cells = 7.22 - 14
    ),
    tolerance = 1e-8
  )
})

test_that(".calculate_ionization_factors: monoprotic base", {
  expect_equal(
    .calculate_ionization_factors(ionization = c("base", 0), pka = c(5, 0)),
    c(
      ion_factor_plasma = 10^(5 - 7.4), # 0.003981072
      ion_factor_cells = 10^(5 - 7.0), # 0.01
      ion_factor_blood_cells = 10^(5 - 7.22) # 0.006025596
    ),
    tolerance = 1e-8
  )
})

test_that(".calculate_ionization_factors: base + acid (zwitterion-like)", {
  # base pKa 5 and acid pKa 7
  expect_equal(
    .calculate_ionization_factors(
      ionization = c("base", "acid"),
      pka = c(5, 7)
    ),
    c(
      ion_factor_plasma = 10^(5 - 7.4) + 10^(7.4 - 7), # 2.515868
      ion_factor_cells = 10^(5 - 7.0) + 10^(7.0 - 7), # 1.01
      ion_factor_blood_cells = 10^(5 - 7.22) + 10^(7.22 - 7) # 1.665613
    ),
    tolerance = 1e-8
  )
})

test_that(".calculate_ionization_factors: a single ionizable group can be in either position", {
  expect_equal(
    .calculate_ionization_factors(
      ionization = c("neutral", "base"),
      pka = c(0, 5)
    ),
    .calculate_ionization_factors(
      ionization = c("base", "neutral"),
      pka = c(5, 0)
    )
  )
  expect_equal(
    .calculate_ionization_factors(
      ionization = c("neutral", "acid"),
      pka = c(0, 7)
    ),
    c(
      ion_factor_plasma = 10^0.4,
      ion_factor_cells = 1,
      ion_factor_blood_cells = 10^0.22
    ),
    tolerance = 1e-8
  )
})
