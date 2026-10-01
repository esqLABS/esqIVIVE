test_that("fit_michaelis_menten_curve fits the expected Km/Vmax from the michaelis_menten_curve.csv fixture", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  mm_curve_path <- system.file(
    "extdata",
    "michaelis_menten_curve.csv",
    package = "ESQivive"
  )
  mm_curve <- read.csv(mm_curve_path)

  fit <- suppressMessages(suppressWarnings(fit_michaelis_menten_curve(
    mm_curve
  )))

  expect_equal(
    fit,
    data.frame(
      parameter = c("km", "vmax"),
      estimate = c(24.819410646102167, 0.148042324167861),
      lower = c(17.526963031269819, 0.130756491082109),
      upper = c(35.521418822256344, 0.170050447145058)
    ),
    tolerance = 1e-6
  )
})
