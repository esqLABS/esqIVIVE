test_that("clearance_IVIVE works", {
  # example for curve
  res <- clearance_IVIVE(
    typeValue = "decay experimentalcurve",
    expData = test_data_cl,
    partitionQSPR = "All Schmitt",
    logLipo = 5,
    ionization = c("acid", 0, 0),
    typeSystem = "hepatocytes",
    FBS_fraction = 0.0,
    pKa = c(3, 0, 0),
    hlcAt_atmm3mol = 1E-6,
    microplateType = 96,
    volMedium_mL = 0.2,
    cCells_Mml = 0.02
  )

  expect_equal(res, 362.3, tolerance = 0.1)

  # example for half life
  expect_equal(
    clearance_IVIVE(
      typeValue = "halfLife",
      expData = 3,
      units = "hours",
      partitionQSPR = "All Schmitt",
      logLipo = 5,
      ionization = c("acid", 0, 0),
      typeSystem = "microsomes",
      FBS_fraction = 0,
      pKa = c(3, 0, 0),
      hlcAt_atmm3mol = 1E-7,
      microplateType = 96,
      volMedium_mL = 0.2,
      cMicro_mgml = 0.2
    ),
    35.9,
    tolerance = 0.1
  )

  # example for clearance values
  expect_equal(
    clearance_IVIVE(
      typeValue = "in vitro clearance parameter",
      expData = 0.3,
      units = "mL/seconds/mg protein",
      partitionQSPR = "All Poulin and Theil",
      logLipo = 2,
      ionization = c("neutral", 0, 0),
      typeSystem = "hepatocytes",
      FBS_fraction = 0,
      pKa = c(0, 0, 0),
      hlcAt_atmm3mol = 0.02,
      microplateType = 96,
      volMedium_mL = 0.2,
      cCells_Mml = 0.2
    ),
    3206.3,
    tolerance = 0.1
  )

  # test clearance with given fu_invitro
  expect_equal(
    clearance_IVIVE(
      typeSystem = "microsomes",
      typeValue = "in vitro clearance parameter",
      expData = 0.3,
      units = "uL/seconds/mg protein",
      partitionQSPR = "All Poulin and Theil",
      logLipo = 2,
      FBS_fraction = 0,
      fu_invitro = 0.1,
      microplateType = 96,
      cMicro_mgml = 0.2
    ),
    9.7,
    tolerance = 0.1
  )
})
