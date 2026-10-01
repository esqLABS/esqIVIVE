test_that(".calculate_ionization_factors: neutral compound", {
  expect_equal(
    .calculate_ionization_factors(ionization = c("neutral", 0), pka = c(0, 0)),
    c(ion_factor_plasma = 0, ion_factor_cells = 0)
  )
})

test_that(".calculate_ionization_factors: monoprotic acid", {
  expect_equal(
    .calculate_ionization_factors(ionization = c("acid", 0), pka = c(14, 0)),
    c(
      ion_factor_plasma = 2.51188643150958e-07,
      ion_factor_cells = 1.65958690743756e-07
    ),
    tolerance = 1e-6
  )
})

test_that(".calculate_ionization_factors: monoprotic base", {
  expect_equal(
    .calculate_ionization_factors(ionization = c("base", 0), pka = c(5, 0)),
    c(
      ion_factor_plasma = 0.00398107170553497,
      ion_factor_cells = 0.00602559586074358
    ),
    tolerance = 1e-6
  )
})

test_that(".calculate_ionization_factors: base + acid (zwitterion-like)", {
  expect_equal(
    .calculate_ionization_factors(
      ionization = c("base", "acid"),
      pka = c(5, 7)
    ),
    c(ion_factor_plasma = 2.51586750321512, ion_factor_cells = 1.6656125032983),
    tolerance = 1e-6
  )
})
