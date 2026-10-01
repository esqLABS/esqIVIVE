test_that("calculate_in_vitro_compartments: no air when the well is full", {
  comp <- calculate_in_vitro_compartments(
    "hepatocytes",
    0,
    96,
    0.392,
    concentration_cells = 0.1
  )

  expect_equal(comp$volume_air, 0)
})

test_that("calculate_in_vitro_compartments gives regression-locked values for hepatocytes, all microplate sizes", {
  expected <- list(
    cell_neutral_lipids = 1.1303e-05,
    cell_neutral_phospholipids = 8.4074e-06,
    cell_acidic_phospholipids = 2.2352e-06,
    cell_proteins = 5.08e-05,
    medium_neutral_lipids = 7.85e-05,
    medium_neutral_phospholipids = 1.5e-05,
    medium_proteins = 0.002
  )

  expected_volume_air <- c(
    "96" = 0.000192,
    "48" = 0.00142,
    "24" = 0.00327,
    "12" = 0.0067
  )
  expected_plastic <- c(
    "96" = 0.606060606060606,
    "48" = 0.363636363636364,
    "24" = 0.257234726688103,
    "12" = 0.181818181818182
  )

  for (mp in c(96, 48, 24, 12)) {
    comp <- calculate_in_vitro_compartments(
      "hepatocytes",
      fbs_fraction = 0.05,
      microplate_type = mp,
      volume_medium = 0.2,
      concentration_cells = 0.1
    )

    expect_equal(
      comp,
      c(
        expected,
        plastic_area_per_volume = expected_plastic[[as.character(mp)]],
        volume_air = expected_volume_air[[as.character(mp)]]
      ),
      tolerance = 1e-6
    )
  }
})

test_that("calculate_in_vitro_compartments gives regression-locked values for microsomes, all microplate sizes", {
  expected <- list(
    cell_neutral_lipids = 0.000261111111111111,
    cell_neutral_phospholipids = 0.000726155555555556,
    cell_acidic_phospholipids = 0.0001594,
    cell_proteins = 0.000740740740740741,
    medium_neutral_lipids = 0,
    medium_neutral_phospholipids = 0,
    medium_proteins = 0,
    plastic_area_per_volume = 0
  )

  expected_volume_air <- c(
    "96" = 0.000192,
    "48" = 0.00142,
    "24" = 0.00327,
    "12" = 0.0067
  )

  for (mp in c(96, 48, 24, 12)) {
    comp <- calculate_in_vitro_compartments(
      "microsomes",
      fbs_fraction = 0,
      microplate_type = mp,
      volume_medium = 0.2,
      concentration_microsomes = 1
    )

    expect_equal(
      comp,
      c(expected, volume_air = expected_volume_air[[as.character(mp)]]),
      tolerance = 1e-6
    )
  }
})

test_that("calculate_in_vitro_compartments rejects invalid inputs", {
  expect_snapshot(
    error = TRUE,
    calculate_in_vitro_compartments(
      "bogus",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.2
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_in_vitro_compartments(
      "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.2
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_in_vitro_compartments(
      "microsomes",
      fbs_fraction = 0,
      microplate_type = 384,
      volume_medium = 0.2,
      concentration_microsomes = 1
    )
  )
})
