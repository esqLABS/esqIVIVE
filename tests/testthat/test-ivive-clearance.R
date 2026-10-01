test_that("ivive_clearance: half_life", {
  expect_equal(
    ivive_clearance(
      value_type = "half_life",
      unit = "hours",
      value = 3,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5
    ),
    2.73522388059702,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance, cells", {
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      species = "human",
      unit = "mL/minutes/millioncells",
      value = 18.27,
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      empirical_correction = FALSE
    ),
    6489.94029850746,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance, microsomes", {
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
    1963343.28358209,
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
    1963343.28358209,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: rate_constant", {
  expect_equal(
    ivive_clearance(
      value_type = "rate_constant",
      unit = "/minutes",
      value = 0.02,
      system = "cells",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5
    ),
    14.2089552238806,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: Wood 2017 empirical correction, 1000-10000 decade", {
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "cells",
      species = "human",
      unit = "mL/minutes/millioncells",
      value = 0.5,
      fu_in_vitro = 0.5,
      concentration_cells = 0.5,
      empirical_correction = TRUE
    ),
    3907.46268656717,
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
      verbose = TRUE
    )
  )
})
