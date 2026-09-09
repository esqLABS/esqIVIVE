test_that("IVIVE_clearance: halfLife", {
  expect_equal(
    IVIVE_clearance(
      typeValue = "halfLife", units = "hours", expData = 3, typeSystem = "hepatocytes",
      fu_invitro = 0.5, cCells_Mml = 0.5
    ),
    2.73522388059702,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance: invitro_clearance_parameter, hepatocytes", {
  expect_equal(
    IVIVE_clearance(
      typeValue = "invitro_clearance_parameter", typeSystem = "hepatocytes", species = "human",
      units = "mL/minutes/millioncells", expData = 18.27, fu_invitro = 0.5, cCells_Mml = 0.5,
      empirical_scalar = "No"
    ),
    6489.94029850746,
    tolerance = 1e-6
  )
})

test_that("IVIVE_clearance: invitro_clearance_parameter, microsomes", {
  expect_equal(
    IVIVE_clearance(
      typeValue = "invitro_clearance_parameter", typeSystem = "microsomes",
      units = "L/minutes/mg protein", expData = 18.27, fu_invitro = 0.5,
      cProtein_mgml = 0.5, volMedium_mL = 0.5, empirical_scalar = "No"
    ),
    1963343.28358209,
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

# KNOWN ISSUE: typeValue="kcat" is accepted by rlang::arg_match() but the
# calculation branch checks for "kcat_min", so it silently falls through to
# ClspePermin <- 0 instead of computing or erroring - not asserted here.
# See plan Known Issues #4.
# KNOWN ISSUE: with empirical_scalar="Yes", the Wood 2017 lookup table has a
# typo ("1000-1000" instead of "1000-10000"), so clearances landing in that
# decade silently become numeric(0) - not asserted here. See plan Known Issues #5.
