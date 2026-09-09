test_that("fit_mm_from_curve fits the expected Km/Vmax from the michaelis_menten_curve.csv fixture", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  mm_curve_path <- system.file("extdata", "michaelis_menten_curve.csv", package = "esqIVIVE")
  mm_curve <- read.csv(mm_curve_path)

  fit <- suppressMessages(suppressWarnings(fit_mm_from_curve(mm_curve)))

  expect_equal(fit["Km_uM", "Mean"], 24.819410646102167, tolerance = 1e-6)
  expect_equal(fit["Vmax_umol_min_mgmicroORcells", "Mean"], 0.148042324167861, tolerance = 1e-6)
})
