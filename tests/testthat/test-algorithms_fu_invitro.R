test_that("calculate_fu_mic_turner", {
  expect_equal(calculate_fu_mic_turner(c("base", 0), c(8, 0), 3, 1), 0.655820505912395, tolerance = 1e-6)
  expect_equal(calculate_fu_mic_turner(c("acid", 0), c(3, 0), 3, 1), 0.897009526377273, tolerance = 1e-6)
  expect_equal(calculate_fu_mic_turner(c("neutral", 0), c(0, 0), 3, 1), 0.574280203654524, tolerance = 1e-6)
})

test_that("calculate_fu_mic_halifax", {
  expect_equal(calculate_fu_mic_halifax(c("base", 0), c(8, 0), 3, 1), 0.795502453049702, tolerance = 1e-6)
  expect_equal(calculate_fu_mic_halifax(c("neutral", 0), c(0, 0), 3, 1), 0.654259613707832, tolerance = 1e-6)
})

test_that("calculate_fu_mic_austin", {
  expect_equal(calculate_fu_mic_austin(c("base", 0), c(8, 0), 3, 1), 0.5689241999032, tolerance = 1e-6)
  expect_equal(calculate_fu_mic_austin(c("neutral", 0), c(0, 0), 3, 1), 0.349395372066127, tolerance = 1e-6)
})

test_that("calculate_fu_hep_austin", {
  expect_equal(
    calculate_fu_hep_austin(ionization = c("base", 0), pKa = c(8, 0), log_lipophilicity = 3, conc_cell_millionml = 0.5),
    0.851936470670237,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_hep_austin(ionization = c("neutral", 0), pKa = c(0, 0), log_lipophilicity = 3, conc_cell_millionml = 0.5),
    0.751683739251381,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_hep_kilford", {
  expect_equal(
    calculate_fu_hep_kilford(ionization = c("base", 0), pKa = c(8, 0), log_lipophilicity = 3, conc_cell_millionml = 0.5),
    0.925640105577411,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_hep_kilford(ionization = c("neutral", 0), pKa = c(0, 0), log_lipophilicity = 3, conc_cell_millionml = 0.5),
    0.858266592080552,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_hep_poulin", {
  expect_equal(
    calculate_fu_hep_poulin(
      ionization = c("base", 0), pKa = c(6, 0), blood_plasma = 1, fraction_unbound = 0.2,
      concentration_cell_neutral_lipids = 0.03, log_lipophilicity = 3
    ),
    0.0334992608857633,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_hep_poulin(
      ionization = c("neutral", 0), pKa = c(0, 0),
      concentration_cell_neutral_lipids = 0.03, log_lipophilicity = 3
    ),
    0.032258064516129,
    tolerance = 1e-6
  )
  # Strong base (ionization[1]=="base" & pKa[1]>7): this branch used to
  # reference an undefined `cCellAPL` and error; it's now a proper parameter.
  expect_equal(
    calculate_fu_hep_poulin(
      ionization = c("base", 0), pKa = c(8, 0), blood_plasma = 1, fraction_unbound = 0.2,
      concentration_cell_neutral_lipids = 0.03, cCellAPL = 0.01, log_lipophilicity = 3
    ),
    0.0203691900058836,
    tolerance = 1e-6
  )
  # cCellAPL is required for the strong-base branch; omitting it errors clearly.
  expect_error(
    calculate_fu_hep_poulin(
      ionization = c("base", 0), pKa = c(8, 0), blood_plasma = 1, fraction_unbound = 0.2,
      concentration_cell_neutral_lipids = 0.03, log_lipophilicity = 3
    )
  )
})
