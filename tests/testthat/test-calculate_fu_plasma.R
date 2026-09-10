test_that("correct_fu_pls_pearce", {
  expect_equal(
    unname(correct_fu_pls_pearce(fraction_unbound = 0.2, log_lipophilicity = 4)),
    0.0133333333333333,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_pls_from_Ks: human", {
  expect_equal(
    unname(calculate_fu_pls_from_Ks(
      partition_albumin = 10^4.48,
      partition_globulin = 10^2.16,
      partition_membrane_lipids = 10^3.51,
      partition_lipids = 100,
      species = "human"
    )),
    0.000798040991411153,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_pls_from_Ks rejects an invalid species", {
  expect_error(
    calculate_fu_pls_from_Ks(
      partition_albumin = 1, partition_globulin = 1, partition_membrane_lipids = 1,
      partition_lipids = 1, species = "bogus"
    )
  )
})
