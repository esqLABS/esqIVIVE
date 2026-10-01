test_that("ivive_clearance: half_life", {
  expect_equal(
    ivive_clearance(
      value_type = "half_life",
      unit = "hours",
      value = 3,
      system = "hepatocytes",
      fu_in_vitro = 0.5,
      concentration_cells = 0.5
    ),
    2.73522388059702,
    tolerance = 1e-6
  )
})

test_that("ivive_clearance: intrinsic_clearance, hepatocytes", {
  expect_equal(
    ivive_clearance(
      value_type = "intrinsic_clearance",
      system = "hepatocytes",
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
      system = "hepatocytes",
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
      system = "hepatocytes",
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
      system = "hepatocytes",
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
      system = "hepatocytes",
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
      system = "hepatocytes"
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
