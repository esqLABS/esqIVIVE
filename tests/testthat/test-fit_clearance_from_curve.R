test_that("fit_clearance_from_curve fits the expected kcat from the clearance.csv fixture", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  kcat <- suppressMessages(suppressWarnings(fit_clearance_from_curve(test_data_cl)))

  expect_equal(unname(kcat["Mean_min-1"]), 0.0187109695257764, tolerance = 1e-6)
  expect_equal(unname(kcat["2.5%_CI"]), 0.0173403013706227, tolerance = 1e-6)
  expect_equal(unname(kcat["95%_CI"]), 0.0202073107873403, tolerance = 1e-6)
})

# KNOWN ISSUE: the R² helper computes `tss` from the time axis (`clear_curve_xy$x`)
# instead of concentration (`$y`), so the reported R² and the "poor fit" warning
# threshold are not meaningful - not asserted here. See plan Known Issues #6.
