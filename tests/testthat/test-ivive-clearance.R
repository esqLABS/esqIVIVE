# Scaling factors in inst/extdata/scaling_factors.csv used below (liver):
#   human: fcell 0.67, MicProtGO 32, CytosProtGO 50, CellsGO 99
#   rat:   fcell 0.72, MicProtGO 44.333, CytosProtGO 73, CellsGO 122.5
# The expected values are written as the scaling formula,
#   clearance = value * unit factor / concentration * (scaling factor / fcell / fu)
# so that the numbers can be followed.

test_that("ivive_clearance: half_life", {
  # 3 hours = 180 min, rate constant = 0.693 / 180 per min, rat hepatocytes
  expect_equal(
    ivive_clearance(
      value_type = "half_life",
      unit = "hours",
      value = 3,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      species = "rat"
    ),
    0.693 / 180 / 0.5 * 122.5 / 0.72 / 0.5,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance, cells", {
  # mL/min/million cells needs no unit conversion and no concentration
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      species = "rat",
      unit = "mL/minutes/millioncells",
      value = 18.27,
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      empirical_correction = FALSE
    ),
    18.27 * 122.5 / 0.72 / 0.5,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance, microsomes", {
  # L/min/mg protein is 1000 times mL/min/mg protein, human microsomes
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "microsomes",
      unit = "L/minutes/mg protein",
      value = 18.27,
      fu_in_vitro = 0.5,
      concentration_microsomes = 0.5,
      volume_medium = 0.5,
      empirical_correction = FALSE
    ),
    18.27 * 1000 * 32 / 0.67 / 0.5,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance per mg protein needs no concentration", {
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "microsomes",
      unit = "L/minutes/mg protein",
      value = 18.27,
      fu_in_vitro = 0.5
    ),
    18.27 * 1000 * 32 / 0.67 / 0.5,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: rate_constant", {
  # rate constant per mL of incubation, divided by the cell concentration
  expect_equal(
    ivive_clearance(
      value_type = "rate_constant",
      unit = "/minutes",
      value = 0.02,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      species = "rat"
    ),
    0.02 / 0.5 * 122.5 / 0.72 / 0.5,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: rate constants per hour and per second are converted to per minute", {
  # 0.02 /min is 1.2 /h and 0.02 / 60 /s
  expected <- 0.02 / 0.5 * 122.5 / 0.72 / 0.5
  expect_equal(
    ivive_clearance(
      value_type = "rate_constant",
      unit = "/hours",
      value = 1.2,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      species = "rat"
    ),
    expected,
    tolerance = 1e-6
  )
  expect_equal(
    ivive_clearance(
      value_type = "rate_constant",
      unit = "/seconds",
      value = 0.02 / 60,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      species = "rat"
    ),
    expected,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: the scaling factors depend on the species", {
  # same incubation, human, rat and dog liver microsomes (fcell, MicProtGO)
  scaled <- function(species) {
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "microsomes",
      unit = "mL/minutes/mg protein",
      value = 1,
      species = species
    )
  }
  expect_equal(scaled("human"), 32 / 0.67, tolerance = 1e-6)
  expect_equal(scaled("rat"), 44.333 / 0.72, tolerance = 1e-6)
  expect_equal(scaled("dog"), 47.525 / 0.72, tolerance = 1e-6)
  expect_equal(scaled("beagle"), 43 / 0.72, tolerance = 1e-6)

  # human hepatocytes
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      unit = "mL/minutes/millioncells",
      value = 1
    ),
    99 / 0.67,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: tissues other than the liver use the average of kidney and gut", {
  # rat microsomes: MicProtGO of the brain is mean(kidney 11, gut 5.8) = 8.4,
  # fcell 0.96
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "microsomes",
      unit = "mL/minutes/mg protein",
      value = 1,
      species = "rat",
      tissue = "brain"
    ),
    mean(c(11, 5.8)) / 0.96,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance per kg uses the organ weight", {
  # human liver weight is 42 g/kg body weight: clearance per kg / weight * 1 / fcell
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "microsomes",
      unit = "mL/minutes/kg",
      value = 100,
      fu_in_vitro = 0.5
    ),
    100 / 42 / 0.67 / 0.5,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: a tissue without a scaling factor gives NA and a warning", {
  # the table has no cells per gram for the human brain (only the liver has)
  expect_snapshot(
    cl <- ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      unit = "mL/minutes/millioncells",
      value = 18.27,
      tissue = "brain"
    )
  )
  expect_true(is.na(cl))

  # nor cytosolic protein for the dog liver, or cells for the dog brain
  expect_warning(
    cl <- ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cytosol",
      unit = "mL/minutes/mg protein",
      value = 1,
      species = "dog"
    ),
    "CytosProtGO"
  )
  expect_true(is.na(cl))
  expect_warning(
    cl <- ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      unit = "mL/minutes/millioncells",
      value = 1,
      species = "dog",
      tissue = "brain"
    ),
    "CellsGO"
  )
  expect_true(is.na(cl))
})

test_that("ivive_clearance: a supported factor of the same tissue gives no warning", {
  # rat and human liver have cells per gram, and human microsomes are
  # available for all tissues
  expect_no_warning(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      unit = "mL/minutes/millioncells",
      value = 1,
      species = "rat"
    )
  )
  expect_no_warning(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "microsomes",
      unit = "mL/minutes/mg protein",
      value = 1,
      tissue = "kidney"
    )
  )
})

test_that("ivive_clearance: the NA of an unsupported scaling factor is kept with the empirical correction", {
  # the empirical correction cannot pick a band for an NA clearance
  expect_warning(
    cl <- ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      unit = "mL/minutes/millioncells",
      value = 1,
      tissue = "brain",
      empirical_correction = TRUE
    ),
    "CellsGO"
  )
  expect_true(is.na(cl))
})

test_that(".apply_wood_correction applies a factor at the band boundaries", {
  # human cells: 0.61, 3.9, 7.1, 22, 1200
  expect_equal(
    .apply_wood_correction(
      10,
      organkgBW = 1,
      system = "cells",
      species = "human"
    ),
    10 * 3.9
  )
  expect_equal(
    .apply_wood_correction(
      100,
      organkgBW = 1,
      system = "cells",
      species = "human"
    ),
    100 * 7.1
  )
  expect_equal(
    .apply_wood_correction(
      1000,
      organkgBW = 1,
      system = "cells",
      species = "human"
    ),
    1000 * 22
  )
  expect_equal(
    .apply_wood_correction(
      10000,
      organkgBW = 1,
      system = "cells",
      species = "human"
    ),
    10000 * 1200
  )
})

test_that("ivive_clearance: Wood 2017 empirical correction, 1000-10000 decade", {
  # rat hepatocytes: the clearance times the liver weight (43 g/kg) is 7316,
  # which is in the 1000-10000 band with a factor of 7.2
  uncorrected <- 0.5 * 122.5 / 0.72 / 0.5
  expect_equal(uncorrected * 43, 7316, tolerance = 1e-4)
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      species = "rat",
      unit = "mL/minutes/millioncells",
      value = 0.5,
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      empirical_correction = TRUE
    ),
    uncorrected * 7.2,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance rejects a unit that does not fit the value type", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "rate_constant",
      unit = "mL/minutes/millioncells",
      value = 18.27,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "intrinsic_clearance",
      unit = "mL/minutes/millioncell",
      value = 18.27,
      system = "cells",
      concentration_cells = 0.5
    )
  )
})

test_that("ivive_clearance rejects a missing concentration", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "minutes",
      value = 3.9,
      system = "microsomes"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "intrinsic_clearance",
      unit = "uL/minutes",
      value = 3.9,
      system = "cells"
    )
  )
})

test_that("ivive_clearance rejects invalid options", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "minutes",
      value = 3.9,
      system = "microsomes",
      concentration_microsomes = 1,
      empirical_correction = "Yes"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "minutes",
      value = 3.9,
      system = "microsomes",
      concentration_microsomes = 1,
      fu_in_vitro = 0
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "minutes",
      value = 3.9,
      system = "microsomes",
      concentration_microsomes = 1,
      species = "dog",
      empirical_correction = TRUE
    )
  )
})

test_that("ivive_clearance rejects an invalid species or tissue", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "hours",
      value = 3,
      system = "microsomes",
      concentration_microsomes = 0.5,
      species = "bogus"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "hours",
      value = 3,
      system = "microsomes",
      concentration_microsomes = 0.5,
      tissue = "bogus"
    )
  )
})

test_that("ivive_clearance: intrinsic_clearance, cytosol uses the cytosolic protein per gram liver", {
  # human liver: CytosProtGO = 50 mg/g, fcell = 0.67
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cytosol",
      unit = "mL/minutes/mg protein",
      value = 0.05,
      fu_in_vitro = 0.8,
      concentration_cytosol = 1
    ),
    0.05 * 50 / 0.67 / 0.8,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: cytosol needs its own concentration", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "minutes",
      value = 3.9,
      system = "cytosol",
      concentration_microsomes = 1
    )
  )
})

test_that("ivive_clearance: the empirical correction is not available for the cytosol", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cytosol",
      unit = "mL/minutes/mg protein",
      value = 0.05,
      concentration_cytosol = 1,
      empirical_correction = TRUE
    )
  )
})

test_that("ivive_clearance rejects the old hepatocytes system", {
  expect_snapshot(
    error = TRUE,
    ivive_clearance(
      value_type = "half_life",
      unit = "hours",
      value = 3,
      system = "hepatocytes",
      concentration_cells = 0.5
    )
  )
})

test_that("ivive_clearance prints only the inputs when verbose", {
  expect_snapshot(
    cl <- ivive_clearance(
      value_type = "intrinsic_clearance",
      value = 18.27,
      unit = "mL/minutes/millioncells",
      system = "cells",
      concentration_cells = 0.5,
      species = "rat",
      verbose = TRUE
    )
  )
})
