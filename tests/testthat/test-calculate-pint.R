test_that("calculate_pint: caco2", {
  expect_equal(
    calculate_pint(method = "caco2", permeability = 2.3E-6),
    6.16744955226714e-06,
    tolerance = 1e-6
  )
})

test_that("calculate_pint: peff", {
  expect_equal(
    calculate_pint(method = "peff", permeability = 2.3E-6),
    1.09553750596963e-09,
    tolerance = 1e-6
  )
})

test_that("calculate_pint rejects an invalid method", {
  expect_snapshot(
    error = TRUE,
    calculate_pint(method = "pampa", permeability = 2.3E-6)
  )
})
