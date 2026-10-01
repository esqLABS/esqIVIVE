test_that("calculate_plasma_partitions: logp method, ionized compound", {
  expect_equal(
    calculate_plasma_partitions(
      method = "logp",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0)
    ),
    list(
      partition_albumin = 0.185104051916421,
      partition_globulin = 0.349000112570181,
      partition_membrane_lipids = 1.000183344634427
    ),
    tolerance = 1e-6
  )
})

test_that("calculate_plasma_partitions: pplfer method", {
  expect_equal(
    calculate_plasma_partitions(
      method = "pplfer",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0),
      lfer_e = 1,
      lfer_b = 0,
      lfer_a = 1.5,
      lfer_s = 0.8,
      lfer_v = 2
    ),
    list(
      partition_albumin = 16.9112310,
      partition_globulin = 4.8520000,
      partition_membrane_lipids = 11324003.6323556
    ),
    tolerance = 1e-6
  )
})

test_that("calculate_plasma_partitions rejects invalid inputs", {
  expect_snapshot(
    error = TRUE,
    calculate_plasma_partitions(
      method = "bogus",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0)
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_plasma_partitions(
      method = "pplfer",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0),
      lfer_e = 1
    )
  )
})

test_that("calculate_plasma_partitions: the ionizable group can be in either position", {
  expect_equal(
    calculate_plasma_partitions(
      method = "logp",
      lipophilicity = 2,
      pka = c(0, 3),
      ionization = c("neutral", "acid")
    ),
    list(
      partition_albumin = 0.185104051916421,
      partition_globulin = 0.349000112570181,
      partition_membrane_lipids = 1.000183344634427
    ),
    tolerance = 1e-6
  )
})
