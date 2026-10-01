test_that("pint_caco2_empir", {
  expect_equal(unname(pint_caco2_empir(Papp_cms = 2.3E-6)), 6.16744955226714e-06, tolerance = 1e-6)
})

test_that("pint_peff_empir", {
  expect_equal(unname(pint_peff_empir(Peff_cms = 2.3E-6)), 1.09553750596963e-09, tolerance = 1e-6)
})
