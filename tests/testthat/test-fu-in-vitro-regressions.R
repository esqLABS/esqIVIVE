test_that(".calculate_fu_mic_turner", {
  expect_equal(
    .calculate_fu_mic_turner(c("base", 0), c(8, 0), 3, 1),
    0.655820505912395,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_mic_turner(c("acid", 0), c(3, 0), 3, 1),
    0.897009526377273,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_mic_turner(c("neutral", 0), c(0, 0), 3, 1),
    0.574280203654524,
    tolerance = 1e-6
  )
})

# Austin et al 2002 and Hallifax and Houston 2006 use logP for strong bases
# (ionization "base" with pKa > 7) and logD at pH 7.4 for acids, weak bases and
# neutral compounds. The expected values below are calculated with the
# regression equations and logD = logP - log10(1 + ion factor), at logP 3 and
# 1 mg/mL microsomal protein.
test_that(".calculate_fu_mic_hallifax", {
  # strong base: logP, the same as a neutral compound
  expect_equal(
    .calculate_fu_mic_hallifax(c("base", 0), c(8, 0), 3, 1),
    0.654259613707832,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_mic_hallifax(c("neutral", 0), c(0, 0), 3, 1),
    0.654259613707832,
    tolerance = 1e-6
  )
  # weak base with pKa 3: almost neutral at pH 7.4, logD = 3 - 1.1e-5
  expect_equal(
    .calculate_fu_mic_hallifax(c("base", 0), c(3, 0), 3, 1),
    0.654264107259219,
    tolerance = 1e-6
  )
  # acid with pKa 3: ionized at pH 7.4, logD = 3 - log10(1 + 10^4.4) = -1.4
  expect_equal(
    .calculate_fu_mic_hallifax(c("acid", 0), c(3, 0), 3, 1),
    0.922994549885457,
    tolerance = 1e-6
  )
})

test_that(".calculate_fu_mic_austin", {
  # strong base: logP, the same as a neutral compound
  expect_equal(
    .calculate_fu_mic_austin(c("base", 0), c(8, 0), 3, 1),
    0.349395372066127,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_mic_austin(c("neutral", 0), c(0, 0), 3, 1),
    0.349395372066127,
    tolerance = 1e-6
  )
  # weak base with pKa 3: almost neutral at pH 7.4
  expect_equal(
    .calculate_fu_mic_austin(c("base", 0), c(3, 0), 3, 1),
    0.349400439815598,
    tolerance = 1e-6
  )
  # acid: logD at pH 7.4. pKa 3 is ionized (logD = -1.4), pKa 6 less so
  expect_equal(
    .calculate_fu_mic_austin(c("acid", 0), c(3, 0), 3, 1),
    0.993643458367834,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_mic_austin(c("acid", 0), c(6, 0), 3, 1),
    0.76948231850321,
    tolerance = 1e-6
  )
})

test_that("the microsomal Austin and Hallifax regressions reproduce the published predictions", {
  # test_fu_microsomes.csv has the logP, class and pKa of 132 measurements with
  # the predictions published for each method (Poulin and Haddad 2011). Calling
  # the regressions with logP reproduces them: logP for strong bases, logD for
  # acids and weak bases, which is the convention of the published methods.
  path <- system.file("extdata", "test_fu_microsomes.csv", package = "ESQivive")
  published <- utils::read.csv(path, check.names = FALSE)
  names(published)[1] <- "Compound"
  pka <- suppressWarnings(as.numeric(published$pKa))
  pka[is.na(pka)] <- 0

  predict <- function(regression) {
    mapply(
      function(class, pka, logp, concentration) {
        regression(c(class, "neutral"), c(pka, 0), logp, concentration)
      },
      published$Class,
      pka,
      published$LogP37C,
      published[["Cp(mg/mL)"]]
    )
  }

  austin <- predict(.calculate_fu_mic_austin)
  hallifax <- predict(.calculate_fu_mic_hallifax)
  difference_austin <- abs(austin - published$Fu_Austin)
  difference_hallifax <- abs(hallifax - published$Fu_HalifaxHouston)

  # one published Austin value is missing
  expect_lt(max(difference_austin, na.rm = TRUE), 0.01)
  expect_lt(max(difference_hallifax), 0.025)

  # by class, since this is where the lipophilicity convention matters
  strong_base <- published$Class == "base" & pka > 7
  acid <- published$Class == "acid"
  expect_lt(max(difference_austin[strong_base], na.rm = TRUE), 0.01)
  expect_lt(max(difference_austin[acid], na.rm = TRUE), 0.01)
  expect_lt(max(difference_hallifax[strong_base]), 0.025)
  expect_lt(max(difference_hallifax[acid]), 0.025)
})

test_that(".calculate_fu_hep_austin", {
  expect_equal(
    .calculate_fu_hep_austin(
      ionization = c("base", 0),
      pka = c(8, 0),
      lipophilicity = 3,
      concentration_cells = 0.5
    ),
    0.851936470670237,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_hep_austin(
      ionization = c("neutral", 0),
      pka = c(0, 0),
      lipophilicity = 3,
      concentration_cells = 0.5
    ),
    0.751683739251381,
    tolerance = 1e-6
  )
})

test_that(".calculate_fu_hep_kilford", {
  expect_equal(
    .calculate_fu_hep_kilford(
      ionization = c("base", 0),
      pka = c(8, 0),
      lipophilicity = 3,
      concentration_cells = 0.5
    ),
    0.925640105577411,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_hep_kilford(
      ionization = c("neutral", 0),
      pka = c(0, 0),
      lipophilicity = 3,
      concentration_cells = 0.5
    ),
    0.858266592080552,
    tolerance = 1e-6
  )
})

test_that(".calculate_fu_hep_poulin", {
  expect_equal(
    .calculate_fu_hep_poulin(
      ionization = c("base", 0),
      pka = c(6, 0),
      lipophilicity = 3,
      concentration_cell_neutral_lipids = 0.03,
      blood_plasma_ratio = 1,
      fu_plasma = 0.2
    ),
    0.0334992608857633,
    tolerance = 1e-6
  )
  expect_equal(
    .calculate_fu_hep_poulin(
      ionization = c("neutral", 0),
      pka = c(0, 0),
      lipophilicity = 3,
      concentration_cell_neutral_lipids = 0.03
    ),
    0.032258064516129,
    tolerance = 1e-6
  )
  # Strong base (ionization[1]=="base" & pka[1]>7): acidic phospholipid binding
  expect_equal(
    .calculate_fu_hep_poulin(
      ionization = c("base", 0),
      pka = c(8, 0),
      lipophilicity = 3,
      concentration_cell_neutral_lipids = 0.03,
      concentration_cell_acidic_phospholipids = 0.01,
      blood_plasma_ratio = 1,
      fu_plasma = 0.2
    ),
    0.0203691900058836,
    tolerance = 1e-6
  )
})
