test_that("fit_depletion_curve fits the expected rate constant from the clearance.csv fixture", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  fit <- suppressMessages(suppressWarnings(fit_depletion_curve(test_data_cl)))

  expect_equal(
    fit,
    data.frame(
      parameter = "rate_constant",
      estimate = 0.0187109695257764,
      lower = 0.0173403013706227,
      upper = 0.0202073107873403
    ),
    tolerance = 1e-6
  )
})
