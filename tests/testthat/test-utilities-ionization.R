test_that("ion_factors: neutral compound", {
  result <- ion_factors(ionization = c("neutral", 0), pKa = c(0, 0))
  expect_equal(unname(result["ion_factor_plasma"]), 0)
  expect_equal(unname(result["ion_factor_cells"]), 0)
})

test_that("ion_factors: monoprotic acid", {
  result <- ion_factors(ionization = c("acid", 0), pKa = c(14, 0))
  expect_equal(unname(result["ion_factor_plasma"]), 2.51188643150958e-07, tolerance = 1e-6)
  expect_equal(unname(result["ion_factor_cells"]), 1.65958690743756e-07, tolerance = 1e-6)
})

test_that("ion_factors: monoprotic base", {
  result <- ion_factors(ionization = c("base", 0), pKa = c(5, 0))
  expect_equal(unname(result["ion_factor_plasma"]), 0.00398107170553497, tolerance = 1e-6)
  expect_equal(unname(result["ion_factor_cells"]), 0.00602559586074358, tolerance = 1e-6)
})

test_that("ion_factors: base + acid (zwitterion-like)", {
  result <- ion_factors(ionization = c("base", "acid"), pKa = c(5, 7))
  expect_equal(unname(result["ion_factor_plasma"]), 2.51586750321512, tolerance = 1e-6)
  expect_equal(unname(result["ion_factor_cells"]), 1.6656125032983, tolerance = 1e-6)
})
