test_that("predict_clearance_ivive works", {
  # example for curve
  res <- predict_clearance_ivive(
    typeValue = "decay experimentalcurve",
    expData = test_data_cl,
    partition_qspr = "All Schmitt",
    log_lipophilicity = 5,
    ionization = c("acid", 0, 0),
    typeSystem = "hepatocytes",
    fetal_bovine_serum_fraction = 0.0,
    pka = c(3, 0, 0),
    henry_law_constant = 1E-6,
    microplate_type = 96,
    volume_medium = 0.2,
    cCells_Mml = 0.02
  )

  expect_equal(res, 362.3, tolerance = 0.1)

  # example for half life
  expect_equal(
    predict_clearance_ivive(
      typeValue = "halfLife",
      expData = 3,
      units = "hours",
      partition_qspr = "All Schmitt",
      log_lipophilicity = 5,
      ionization = c("acid", 0, 0),
      typeSystem = "microsomes",
      fetal_bovine_serum_fraction = 0,
      pka = c(3, 0, 0),
      henry_law_constant = 1E-7,
      microplate_type = 96,
      volume_medium = 0.2,
      cMicro_mgml = 0.2
    ),
    35.9,
    tolerance = 0.1
  )

  # example for clearance values
  expect_equal(
    predict_clearance_ivive(
      typeValue = "in vitro clearance parameter",
      expData = 0.3,
      units = "mL/seconds/mg protein",
      partition_qspr = "All Poulin and Theil",
      log_lipophilicity = 2,
      ionization = c("neutral", 0, 0),
      typeSystem = "hepatocytes",
      fetal_bovine_serum_fraction = 0,
      pka = c(0, 0, 0),
      henry_law_constant = 0.02,
      microplate_type = 96,
      volume_medium = 0.2,
      cCells_Mml = 0.2
    ),
    3206.3,
    tolerance = 0.1
  )

  # test clearance with given fu_invitro
  expect_equal(
    predict_clearance_ivive(
      typeSystem = "microsomes",
      typeValue = "in vitro clearance parameter",
      expData = 0.3,
      units = "uL/seconds/mg protein",
      partition_qspr = "All Poulin and Theil",
      log_lipophilicity = 2,
      fetal_bovine_serum_fraction = 0,
      fu_invitro = 0.1,
      microplate_type = 96,
      cMicro_mgml = 0.2
    ),
    9.7,
    tolerance = 0.1
  )
})
