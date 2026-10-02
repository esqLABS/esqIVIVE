# Scaling factors in inst/extdata/scaling_factors.csv used below (liver):
#   human: fcell 0.67, MicProtGO 32, CellsGO 99
#   rat:   fcell 0.72, MicProtGO 44.333, CellsGO 122.5
# vmax in vivo = vmax * scaling factor / fcell * 1000 g/L, km_unbound = km * fu

test_that("ivive_michaelis_menten: cells", {
  expect_equal(
    ivive_michaelis_menten(
      system = "cells",
      vmax = 2,
      km = 1,
      tissue = "liver",
      species = "rat",
      relative_expression_factor = 1
    ),
    list(vmax = 2 * 122.5 / 0.72 * 1000, km_unbound = 1),
    tolerance = 1e-6
  )
})

test_that("ivive_michaelis_menten: microsomes", {
  expect_equal(
    ivive_michaelis_menten(
      system = "microsomes",
      fu_in_vitro = 0.2,
      vmax = 2,
      km = 1
    ),
    list(vmax = 2 * 32 / 0.67 * 1000, km_unbound = 0.2),
    tolerance = 1e-6
  )
})

test_that("ivive_michaelis_menten: the relative expression factor scales vmax", {
  expect_equal(
    ivive_michaelis_menten(
      system = "microsomes",
      vmax = 2,
      km = 1,
      species = "dog",
      relative_expression_factor = 0.5
    )$vmax,
    2 * 47.525 * 0.5 / 0.72 * 1000,
    tolerance = 1e-6
  )
})

test_that("ivive_michaelis_menten: a tissue without a scaling factor gives NA vmax and a warning", {
  # the table has no cells per gram for the human brain (only the liver has)
  expect_snapshot(
    mm <- ivive_michaelis_menten(
      system = "cells",
      vmax = 2,
      km = 1,
      tissue = "brain"
    )
  )
  expect_true(is.na(mm$vmax))
  # the unbound Km does not need a scaling factor
  expect_equal(mm$km_unbound, 1)

  expect_warning(
    mm <- ivive_michaelis_menten(
      system = "cells",
      vmax = 2,
      km = 1,
      species = "dog",
      tissue = "lung"
    ),
    "CellsGO"
  )
  expect_true(is.na(mm$vmax))
})

test_that("ivive_michaelis_menten: human hepatocytes", {
  expect_equal(
    ivive_michaelis_menten(system = "cells", vmax = 2, km = 1)$vmax,
    2 * 99 / 0.67 * 1000,
    tolerance = 1e-6
  )
})

test_that("ivive_michaelis_menten: no warning when the scaling factor is available", {
  expect_no_warning(
    ivive_michaelis_menten(system = "cells", vmax = 2, km = 1, species = "rat")
  )
  expect_no_warning(
    ivive_michaelis_menten(
      system = "microsomes",
      vmax = 2,
      km = 1,
      species = "beagle",
      tissue = "gut"
    )
  )
})

test_that("ivive_michaelis_menten rejects invalid inputs", {
  expect_snapshot(
    error = TRUE,
    ivive_michaelis_menten(
      system = "cells",
      vmax = 2,
      km = 1,
      species = "bogus"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_michaelis_menten(
      system = "cells",
      vmax = 2,
      km = 1,
      tissue = "bogus"
    )
  )
  expect_snapshot(
    error = TRUE,
    ivive_michaelis_menten(
      system = "cells",
      vmax = 2,
      km = 1,
      fu_in_vitro = 1.2
    )
  )
})
