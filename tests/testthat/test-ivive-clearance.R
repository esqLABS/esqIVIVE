test_that("IVIVE_clearance: halfLife", {
  expect_equal(
    unname(IVIVE_clearance(
      typeValue = "halfLife", units = "hours", expData = 3, typeSystem = "hepatocytes",
      fu_invitro = 0.5, cCells_Mml = 0.5
    )),
    2.73522388059702,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance: invitro_clearance_parameter, hepatocytes", {
  expect_equal(
    unname(IVIVE_clearance(
      typeValue = "invitro_clearance_parameter", typeSystem = "hepatocytes", species = "human",
      units = "mL/minutes/millioncells", expData = 18.27, fu_invitro = 0.5, cCells_Mml = 0.5,
      empirical_scalar = "No"
    )),
    6489.94029850746,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance: invitro_clearance_parameter, microsomes", {
  expect_equal(
    unname(IVIVE_clearance(
      typeValue = "invitro_clearance_parameter", typeSystem = "microsomes",
      units = "L/minutes/mg protein", expData = 18.27, fu_invitro = 0.5,
      cProtein_mgml = 0.5, volMedium_mL = 0.5, empirical_scalar = "No"
    )),
    1963343.28358209,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance: kcat, correct rate-constant units", {
  expect_equal(
    unname(IVIVE_clearance(
      typeValue = "kcat", units = "/minutes", expData = 0.02, typeSystem = "hepatocytes",
      fu_invitro = 0.5, cCells_Mml = 0.5
    )),
    14.2089552238806,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance: kcat warns when given non-rate-constant units", {
  expect_warning(
    IVIVE_clearance(
      typeValue = "kcat", units = "mL/minutes/millioncells", expData = 18.27, typeSystem = "hepatocytes",
      fu_invitro = 0.5, cCells_Mml = 0.5
    )
  )
})

test_that("IVIVE_clearance: Wood 2017 empirical scalar, 1000-10000 decade", {
  expect_equal(
    unname(IVIVE_clearance(
      typeValue = "invitro_clearance_parameter", typeSystem = "hepatocytes", species = "human",
      units = "mL/minutes/millioncells", expData = 0.5, fu_invitro = 0.5, cCells_Mml = 0.5,
      empirical_scalar = "Yes"
    )),
    3907.46268656717,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance rejects an invalid species or tissue", {
  expect_error(
    IVIVE_clearance(
      typeValue = "halfLife", units = "hours", expData = 3, typeSystem = "microsomes",
      cProtein_mgml = 0.5, species = "bogus"
    )
  )
  expect_error(
    IVIVE_clearance(
      typeValue = "halfLife", units = "hours", expData = 3, typeSystem = "microsomes",
      cProtein_mgml = 0.5, tissue = "bogus"
    )
  )
})
