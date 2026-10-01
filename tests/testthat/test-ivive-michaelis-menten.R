test_that("ivive_michaelis_menten: hepatocytes", {
  expect_equal(
    ivive_michaelis_menten(
      system = "hepatocytes",
      vmax = 2,
      km = 1,
      tissue = "liver",
      species = "human",
      relative_expression_factor = 1
    ),
    list(vmax = 355223.880597015, km_unbound = 1),
    tolerance = 1e-6
  )
})

test_that("ivive_michaelis_menten: microsomes", {
  expect_equal(
    ivive_michaelis_menten(
      system = "microsomes",
      fu_in_vitro = 0.2,
      vmax = 2,
      km = 1
    ),
    list(vmax = 107462.686567164, km_unbound = 0.2),
    tolerance = 1e-6
  )
})

test_that("ivive_michaelis_menten rejects invalid inputs", {
  expect_snapshot(
    error = TRUE,
    ivive_michaelis_menten(
      system = "hepatocytes",
      vmax = 2,
      km = 1,
      species = "bogus"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_michaelis_menten(
      system = "hepatocytes",
      vmax = 2,
      km = 1,
      tissue = "bogus"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_michaelis_menten(
      system = "hepatocytes",
      vmax = 2,
      km = 1,
      fu_in_vitro = 1.2
    )
  )
})
