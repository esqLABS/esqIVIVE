test_that("correct_fu_plasma_pearce", {
  expect_equal(
    correct_fu_plasma_pearce(fu_plasma = 0.2, lipophilicity = 4),
    0.0133333333333333,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_plasma: human", {
  expect_equal(
    calculate_fu_plasma(
      partition_albumin = 10^4.48,
      partition_globulin = 10^2.16,
      partition_membrane_lipids = 10^3.51,
      partition_neutral_lipids = 100,
      species = "human"
    ),
    0.000798040991411153,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_plasma rejects an invalid species", {
  expect_snapshot(
    error = TRUE,
    calculate_fu_plasma(
      partition_albumin = 1,
      partition_globulin = 1,
      partition_membrane_lipids = 1,
      partition_neutral_lipids = 1,
      species = "bogus"
    )
  )
})
