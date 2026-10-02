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

test_that("fit_depletion_curve leaves out rows with a missing concentration", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  incomplete <- rbind(
    test_data_cl,
    data.frame(Time_min = 30, Concentration_uM = NA)
  )
  names(incomplete) <- names(test_data_cl)

  expect_equal(
    suppressMessages(suppressWarnings(fit_depletion_curve(incomplete))),
    data.frame(
      parameter = "rate_constant",
      estimate = 0.0187109695257764,
      lower = 0.0173403013706227,
      upper = 0.0202073107873403
    ),
    tolerance = 1e-6
  )
})
