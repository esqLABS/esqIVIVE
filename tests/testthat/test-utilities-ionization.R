# Reference values follow from the Henderson-Hasselbalch relations, with
# plasma pH 7.4, intracellular pH 7.0 and blood cell pH 7.22:
#   acid: ion factor = 10^(pH - pKa),  base: ion factor = 10^(pKa - pH)
# A base + acid compound is the sum of both contributions.
# An ion factor is the ratio of the ionized to the neutral form: 0 means fully
# neutral, 1 means half ionized, and large values mean mostly ionized.

test_that(".calculate_ionization_factors: neutral compound", {
  # Tests that a compound without an ionizable group has no ionized form, so
  # all three factors are exactly 0 whatever the pKa.
  expect_equal(
    .calculate_ionization_factors(ionization = c("neutral", 0), pka = c(0, 0)),
    c(ion_factor_plasma = 0, ion_factor_cells = 0, ion_factor_blood_cells = 0)
  )
})

test_that(".calculate_ionization_factors: monoprotic acid", {
  # Tests the acid formula at the three pH values (plasma, cells, blood
  # cells). A pKa of 7 is used so the factors are close to 1 and the
  # comparison is meaningful: in cells (pH 7.0) the acid is exactly half
  # ionized, so its factor is 1.
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
  # Tests an acid that is hardly ionized (pKa 14), where the factors are
  # about 1e-7. They are compared on the log10 scale because testthat falls
  # back to an absolute tolerance for values this small, which would let
  # almost any number pass.
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
  # Tests the base formula, which is the mirror image of the acid one: a base
  # is more ionized at lower pH, so the factor is 10^(pKa - pH). With a pKa of
  # 5 the factor is largest in the most acidic compartment, the cells.
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
  # Tests a compound with one basic group (pKa 5) and one acidic group
  # (pKa 7): the factor at each pH is the sum of the base and the acid
  # contributions.
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
  # Tests that a single ionizable group gives the same factors in the second
  # slot (c("neutral", "base"), pKa c(0, 5)) as in the first (c("base",
  # "neutral"), pKa c(5, 0)), for a base and for an acid.
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
