test_that("generating an in vitro compartment works", {
  comp <- in_vitro_compartments(
    "hepatocytes",
    0,
    96,
    0.392,
    cCells_Mml = 0.1
  )

  expect_equal(comp$volAir_L, 0)
})

test_that("in_vitro_compartments gives regression-locked values for hepatocytes, all microplate sizes", {
  expected <- list(
    cCellNL_vvmedium = 1.1303e-05,
    cCellNPL_vvmedium = 8.4074e-06,
    cCellAPL_vvmedium = 2.2352e-06,
    cCellPro_vvmedium = 5.08e-05,
    cMediumNL_vvmedium = 7.85e-05,
    cMediumNPL_vvmedium = 1.5e-05,
    cMediumPro_vvmedium = 0.002
  )

  expected_volAir <- c("96" = 0.000192, "48" = 0.00142, "24" = 0.00327, "12" = 0.0067)
  expected_saPlastic <- c(
    "96" = 0.606060606060606,
    "48" = 0.363636363636364,
    "24" = 0.257234726688103,
    "12" = 0.181818181818182
  )

  for (mp in c(96, 48, 24, 12)) {
    comp <- in_vitro_compartments(
      "hepatocytes",
      FBS_fraction = 0.05,
      microplateType = mp,
      volMedium_mL = 0.2,
      cCells_Mml = 0.1
    )

    for (nm in names(expected)) {
      expect_equal(comp[[nm]], expected[[nm]], tolerance = 1e-6)
    }
    expect_equal(comp$volAir_L, expected_volAir[[as.character(mp)]], tolerance = 1e-6)
    expect_equal(
      comp$saPlasticVolMedium_m2L,
      expected_saPlastic[[as.character(mp)]],
      tolerance = 1e-6
    )
  }
})

test_that("in_vitro_compartments gives regression-locked values for microsomes, all microplate sizes", {
  expected <- list(
    cCellNL_vvmedium = 0.000261111111111111,
    cCellNPL_vvmedium = 0.000726155555555556,
    cCellAPL_vvmedium = 0.0001594,
    cCellPro_vvmedium = 0.000740740740740741,
    cMediumNL_vvmedium = 0,
    cMediumNPL_vvmedium = 0,
    cMediumPro_vvmedium = 0,
    saPlasticVolMedium_m2L = 0
  )

  expected_volAir <- c("96" = 0.000192, "48" = 0.00142, "24" = 0.00327, "12" = 0.0067)

  for (mp in c(96, 48, 24, 12)) {
    comp <- in_vitro_compartments(
      "microsomes",
      FBS_fraction = 0,
      microplateType = mp,
      volMedium_mL = 0.2,
      cMicro_mgml = 1
    )

    for (nm in names(expected)) {
      expect_equal(comp[[nm]], expected[[nm]], tolerance = 1e-6)
    }
    expect_equal(comp$volAir_L, expected_volAir[[as.character(mp)]], tolerance = 1e-6)
  }
})

test_that("in_vitro_compartments rejects an invalid typeSystem", {
  expect_error(
    in_vitro_compartments("bogus", FBS_fraction = 0, microplateType = 96, volMedium_mL = 0.2)
  )
})
