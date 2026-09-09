test_that("calculate_fu_in_vitro: All PK-Sim Standard", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "All PK-Sim Standard", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2
    ),
    0.468271768293341,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: PK-Sim Standard + fu", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "PK-Sim Standard + fu", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "hepatocytes", FBS_fraction = 0.05, microplate_type = 96, fraction_unbound = 0.2,
      volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2
    ),
    0.428171630942448,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: Poulin and Theil + fu", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin and Theil + fu", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96, fraction_unbound = 0.01, blood_plasma_ratio = 2,
      volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2
    ),
    0.783304715148699,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: All Schmitt", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "All Schmitt", log_lipophilicity = 0.42, ionization = c("acid", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96, fraction_unbound = 0.2, blood_plasma_ratio = 1,
      volume_medium = 0.22, pka = c(6, 0), concentration_microsomes = 2
    ),
    0.99122141082607,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: Schmitt + fu", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Schmitt + fu", log_lipophilicity = 0.42, ionization = c("acid", 0),
      type_system = "microsomes", FBS_fraction = 0.05, microplate_type = 96, fraction_unbound = 0.2, blood_plasma_ratio = 1,
      volume_medium = 0.22, pka = c(6, 0), concentration_microsomes = 2
    ),
    0.827228158380981,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: literature methods, microsomes", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Austin", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(3, 0), concentration_microsomes = 1
    ),
    0.349395372066127,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Halifax", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(3, 0), concentration_microsomes = 1
    ),
    0.654259613707832,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Turner", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(3, 0), concentration_microsomes = 1
    ),
    0.897009526377273,
    tolerance = 1e-6
  )
  # Poulin with a neutral compound, to avoid the known undefined-cCellAPL issue
  # in the strong-base (pKa>7) branch of calculate_fu_hep_poulin() - see Known Issues.
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin", log_lipophilicity = 3, ionization = c("neutral", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), concentration_microsomes = 1
    ),
    0.503203730416988,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: literature methods, hepatocytes", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Austin", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(3, 0), concentration_cells = 1
    ),
    0.602158093174717,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Kilford", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(3, 0), concentration_cells = 1
    ),
    0.751722412716774,
    tolerance = 1e-6
  )
  # Poulin with a neutral compound, to avoid the known undefined-cCellAPL issue - see Known Issues.
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin", log_lipophilicity = 3, ionization = c("neutral", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), concentration_cells = 1
    ),
    0.835349309667331,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: All_literature (neutral compound, avoids the Poulin base>7 issue)", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "All_literature", log_lipophilicity = 3, ionization = c("neutral", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), concentration_microsomes = 1
    ),
    0.520284729961368,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "All_literature", log_lipophilicity = 3, ionization = c("neutral", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), concentration_cells = 1
    ),
    0.72974327185294,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro rejects invalid partition_qspr and type_system", {
  expect_error(
    calculate_fu_in_vitro(
      partition_qspr = "bogus", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(6, 0), concentration_cells = 2
    )
  )
  expect_error(
    calculate_fu_in_vitro(
      partition_qspr = "All PK-Sim Standard", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "bogus", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(6, 0), concentration_cells = 2
    )
  )
})

# KNOWN ISSUE: "Rodgers & Rowland + fu" uses an undefined variable `X` and always
# errors - not asserted here. See plan Known Issues #1.
# KNOWN ISSUE: "All Schmitt"/"Schmitt + fu" with a two-group ionization (e.g.
# c("acid","base")) hits a copy-paste bug in calculate_ionization_schmitt()'s F3
# guard - not asserted here. See plan Known Issues #7.
# KNOWN ISSUE: the Poulin literature branch with a strong base (pKa[1] > 7)
# references an undefined `cCellAPL` inside calculate_fu_hep_poulin() and errors -
# not asserted here. See plan Known Issues #2.
