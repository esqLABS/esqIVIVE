test_that("generating an in vitro compartment works", {
  comp <- calculate_in_vitro_compartments(
    "hepatocytes",
    0,
    96,
    0.392,
    cCells_Mml = 0.1
  )

  expect_equal(comp$volAir_L, 0)
})
