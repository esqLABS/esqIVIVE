test_that("clearance_IVIVE works", {
  # example for curve
  res <- clearance_IVIVE(
    typeValue = "decay experimentalcurve",
    expData = test_data_cl,
    partitionQSPR = "All Schmitt",
    logLipo = 5,
    ionization = c("acid", 0, 0),
    typeSystem = "hepatocytes",
    FBS = 0.0, pKa = c(3, 0, 0),
    hlcAt = 1E-6,
    microplateType = 96,
    volMedium = 0.2,
    cCells = 0.02
  )

  expect_equal(res, 283.6378, tolerance = 0.0001)

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
      FBS = 0, pKa = c(3, 0, 0),
      hlcAt = 1E-7,
      microplateType = 96,
      volMedium = 0.2,
      cMicro = 0.2
    ),
    27.8836,
    tolerance = 0.0001
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
      FBS = 0, pKa = c(0, 0, 0),
      hlcAt = 0.02,
      microplateType = 96,
      volMedium = 0.2,
      cCells = 0.2
    ),
    3745.149,
    tolerance = 0.0001
  )

  # test clearance with given fu_hep
  expect_equal(
    clearance_IVIVE(
      typeSystem = "microsomes",
      typeValue = "in vitro clearance parameter",
      expData = 0.3,
      units = "uL/seconds/mg protein",
      partitionQSPR = "All Poulin and Theil",
      logLipo = 2,
      FBS = 0, fu_hep = 0.1,
      microplateType = 96,
      cMicro = 0.2
    ),
    10.746,
    tolerance = 0.0001
  )
})
