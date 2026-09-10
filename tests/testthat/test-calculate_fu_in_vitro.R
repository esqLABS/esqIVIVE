test_that("calculate_fu_in_vitro: All PK-Sim Standard", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "All PK-Sim Standard", log_lipophilicity = 3, ionization = c("acid", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2
    ),
    0.563008780871692,
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
    0.506029024084582,
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
    0.783309951145294,
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
    0.992174983463575,
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
    0.827892197909482,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: All Schmitt, zwitterion (two ionizable groups)", {
  # calculate_ionization_schmitt() used to return NA for any compound with both
  # ionization groups set (F3 copy-paste bug); the unused three-group scaffolding
  # was removed and this now computes a real number.
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "All Schmitt", log_lipophilicity = 2, ionization = c("acid", "base"),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(4, 9), concentration_microsomes = 1
    ),
    0.876012848119017,
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
  # Poulin, neutral compound (else branch of calculate_fu_hep_poulin's ionization check)
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin", log_lipophilicity = 3, ionization = c("neutral", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), concentration_microsomes = 1
    ),
    0.503203730416988,
    tolerance = 1e-6
  )
  # Poulin, strong base (branch that used to error on an undefined cCellAPL - now fixed)
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(8, 0), concentration_microsomes = 1,
      blood_plasma_ratio = 1, fraction_unbound = 0.2
    ),
    0.535098540744975,
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
  # Poulin, neutral compound (else branch of calculate_fu_hep_poulin's ionization check)
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin", log_lipophilicity = 3, ionization = c("neutral", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), concentration_cells = 1
    ),
    0.835349309667331,
    tolerance = 1e-6
  )
  # Poulin, strong base (branch that used to error on an undefined cCellAPL - now fixed)
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Poulin", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(8, 0), concentration_cells = 1,
      blood_plasma_ratio = 1, fraction_unbound = 0.2
    ),
    0.88213946637191,
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

test_that("calculate_fu_in_vitro: Rodgers & Rowland + fu (strong base only)", {
  expect_equal(
    calculate_fu_in_vitro(
      partition_qspr = "Rodgers & Rowland + fu", log_lipophilicity = 3, ionization = c("base", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(8, 0), henry_law_constant = 1E-6,
      fraction_unbound = 0.2, blood_plasma_ratio = 1, concentration_cells = 2
    ),
    0.94721416433249,
    tolerance = 1e-6
  )

  # Only strong bases (ionization[1]=="base" & pka[1]>7) are supported: PK-Sim
  # uses a different, protein-binding-based pathway (Ka_PR) for acids, neutrals
  # and weak bases that this branch does not implement.
  expect_error(
    calculate_fu_in_vitro(
      partition_qspr = "Rodgers & Rowland + fu", log_lipophilicity = 2, ionization = c("neutral", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(0, 0), henry_law_constant = 1E-6,
      fraction_unbound = 0.3, blood_plasma_ratio = 1, concentration_cells = 2
    )
  )
  expect_error(
    calculate_fu_in_vitro(
      partition_qspr = "Rodgers & Rowland + fu", log_lipophilicity = 1, ionization = c("acid", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(4, 0), henry_law_constant = 1E-6,
      fraction_unbound = 0.5, blood_plasma_ratio = 0.8, concentration_cells = 2
    )
  )
  expect_error(
    calculate_fu_in_vitro(
      partition_qspr = "Rodgers & Rowland + fu", log_lipophilicity = 2, ionization = c("base", 0),
      type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
      volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6,
      fraction_unbound = 0.3, blood_plasma_ratio = 1, concentration_cells = 2
    )
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
