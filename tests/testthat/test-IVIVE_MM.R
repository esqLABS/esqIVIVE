test_that("IVIVE_MM: hepatocytes", {
  result <- IVIVE_MM(typeSystem = "hepatocytes", vmax = 2, km_micromolar = 1, tissue = "liver", species = "human", REF = 1)

  expect_equal(result$vmax_umol_minL, 355223.880597015, tolerance = 1e-6)
  expect_equal(result$Km_unb_uM, 1, tolerance = 1e-6)
})

test_that("IVIVE_MM: microsomes", {
  result <- IVIVE_MM(typeSystem = "microsomes", fu_invitro = 0.2, vmax = 2, km_micromolar = 1)

  expect_equal(result$vmax_umol_minL, 107462.686567164, tolerance = 1e-6)
  expect_equal(result$Km_unb_uM, 0.2, tolerance = 1e-6)
})

test_that("IVIVE_MM rejects an invalid species or tissue", {
  expect_error(IVIVE_MM(typeSystem = "hepatocytes", vmax = 2, km_micromolar = 1, species = "bogus"))
  expect_error(IVIVE_MM(typeSystem = "hepatocytes", vmax = 2, km_micromolar = 1, tissue = "bogus"))
})
