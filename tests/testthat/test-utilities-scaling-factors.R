scaling_factors <- utils::read.csv(
  system.file("extdata", "scaling_factors.csv", package = "ESQivive")
)

test_that("the scaling factors of tissues other than the liver are the average of kidney and gut", {
  # documented in ivive_clearance() and ivive_michaelis_menten(): the factors
  # per gram tissue of the other tissues are the mean of kidney and gut, and
  # are NA if one of the two is NA. The literature values are rounded, so a
  # tolerance of 0.2 is used.
  per_gram <- c("MicProtGO", "CytosProtGO", "CellsGO")
  for (species in unique(scaling_factors$species)) {
    rows <- scaling_factors[scaling_factors$species == species, ]
    row <- function(organ) rows[rows$organ == organ, per_gram]
    expected <- (unlist(row("kidney")) + unlist(row("gut"))) / 2

    others <- rows[!rows$organ %in% c("liver", "kidney", "gut"), ]
    for (organ in others$organ) {
      actual <- unlist(row(organ))
      expect_equal(
        is.na(actual),
        is.na(expected),
        info = paste(species, organ, "NA pattern")
      )
      expect_equal(
        actual[!is.na(actual)],
        expected[!is.na(expected)],
        tolerance = 0.2 / 30,
        info = paste(species, organ)
      )
    }
  }
})

test_that("the scaling factor table has the species and tissues of the IVIVE functions", {
  expect_setequal(
    unique(scaling_factors$species),
    c("human", "rat", "dog", "beagle")
  )
  expect_setequal(
    unique(scaling_factors$organ),
    c("brain", "gonads", "heart", "kidney", "gut", "liver", "lung")
  )
  # every species has every tissue once, with fcell and the organ weight
  expect_equal(nrow(scaling_factors), 4 * 7)
  expect_false(anyNA(scaling_factors[c("fcell", "weightorgankgBW")]))
})

test_that(".warn_unsupported_scaling_factors only warns for NA factors", {
  row <- data.frame(fcell = 0.5, MicProtGO = NA_real_, CellsGO = 10)
  expect_no_warning(
    .warn_unsupported_scaling_factors(row, c("fcell", "CellsGO"), "human", "liver")
  )
  expect_snapshot(
    .warn_unsupported_scaling_factors(row, "MicProtGO", "human", "brain")
  )
})
