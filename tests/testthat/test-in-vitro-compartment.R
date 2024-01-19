test_that("generating an in vitro compartment works", {
  comp <- getInVitroCompartment(1, 1, 96, 0.392)

  expect_equal(comp$volAir_L, 0)
})
