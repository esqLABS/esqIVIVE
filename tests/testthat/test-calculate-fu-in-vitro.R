test_that("calculate_fu_in_vitro: pksim_standard", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "pksim_standard",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(6, 0),
      henry_law_constant = 1E-6,
      concentration_cells = 2
    ),
    0.563006547940565,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: pksim_standard_fu", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "pksim_standard_fu",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "cells",
      fbs_fraction = 0.05,
      microplate_type = 96,
      fu_plasma = 0.2,
      volume_medium = 0.22,
      pka = c(6, 0),
      henry_law_constant = 1E-6,
      concentration_cells = 2
    ),
    0.506027220252861,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: poulin_theil_fu", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "poulin_theil_fu",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      fu_plasma = 0.01,
      blood_plasma_ratio = 2,
      volume_medium = 0.22,
      pka = c(6, 0),
      henry_law_constant = 1E-6,
      concentration_cells = 2
    ),
    0.783305628878156,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: schmitt", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "schmitt",
      lipophilicity = 0.42,
      ionization = c("acid", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      fu_plasma = 0.2,
      blood_plasma_ratio = 1,
      volume_medium = 0.22,
      pka = c(6, 0),
      concentration_microsomes = 2
    ),
    0.992174983463575,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: schmitt_fu", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "schmitt_fu",
      lipophilicity = 0.42,
      ionization = c("acid", 0),
      system = "microsomes",
      fbs_fraction = 0.05,
      microplate_type = 96,
      fu_plasma = 0.2,
      blood_plasma_ratio = 1,
      volume_medium = 0.22,
      pka = c(6, 0),
      concentration_microsomes = 2
    ),
    0.827892197909482,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: schmitt, zwitterion (two ionizable groups)", {
  # .calculate_ionization_schmitt() used to return NA for any compound with both
  # ionization groups set (F3 copy-paste bug); the unused three-group scaffolding
  # was removed and this now computes a real number.
  expect_equal(
    calculate_fu_in_vitro(
      method = "schmitt",
      lipophilicity = 2,
      ionization = c("acid", "base"),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(4, 9),
      concentration_microsomes = 1
    ),
    0.876012848119017,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: literature methods, microsomes", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "austin",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(3, 0),
      concentration_microsomes = 1
    ),
    0.349395372066127,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      method = "hallifax",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(3, 0),
      concentration_microsomes = 1
    ),
    0.654259613707832,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      method = "turner",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(3, 0),
      concentration_microsomes = 1
    ),
    0.897009526377273,
    tolerance = 1e-6
  )
  # Poulin, neutral compound (else branch of the Poulin ionization check)
  expect_equal(
    calculate_fu_in_vitro(
      method = "poulin",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(0, 0),
      concentration_microsomes = 1
    ),
    0.503203730416988,
    tolerance = 1e-6
  )
  # Poulin, strong base
  expect_equal(
    calculate_fu_in_vitro(
      method = "poulin",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(8, 0),
      concentration_microsomes = 1,
      blood_plasma_ratio = 1,
      fu_plasma = 0.2
    ),
    0.535098540744975,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: literature methods, cells", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "austin",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(3, 0),
      concentration_cells = 1
    ),
    0.602158093174717,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      method = "kilford",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(3, 0),
      concentration_cells = 1
    ),
    0.751722412716774,
    tolerance = 1e-6
  )
  # Poulin, neutral compound (else branch of the Poulin ionization check)
  expect_equal(
    calculate_fu_in_vitro(
      method = "poulin",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(0, 0),
      concentration_cells = 1
    ),
    0.835349309667331,
    tolerance = 1e-6
  )
  # Poulin, strong base
  expect_equal(
    calculate_fu_in_vitro(
      method = "poulin",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(8, 0),
      concentration_cells = 1,
      blood_plasma_ratio = 1,
      fu_plasma = 0.2
    ),
    0.88213946637191,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: all_literature, neutral compound", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "all_literature",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      system = "microsomes",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(0, 0),
      concentration_microsomes = 1
    ),
    0.520284729961368,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      method = "all_literature",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(0, 0),
      concentration_cells = 1
    ),
    0.72974327185294,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro: rodgers_rowland_fu, strong bases only", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "rodgers_rowland_fu",
      lipophilicity = 3,
      ionization = c("base", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(8, 0),
      henry_law_constant = 1E-6,
      fu_plasma = 0.2,
      blood_plasma_ratio = 1,
      concentration_cells = 2
    ),
    0.947223702721102,
    tolerance = 1e-6
  )

  # Only strong bases (ionization[1]=="base" & pka[1]>7) are supported: PK-Sim
  # uses a different, protein-binding-based pathway (Ka_PR) for acids, neutrals
  # and weak bases that this branch does not implement.
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "rodgers_rowland_fu",
      lipophilicity = 2,
      ionization = c("neutral", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(0, 0),
      henry_law_constant = 1E-6,
      fu_plasma = 0.3,
      blood_plasma_ratio = 1,
      concentration_cells = 2
    )
  )
  expect_error(
    calculate_fu_in_vitro(
      method = "rodgers_rowland_fu",
      lipophilicity = 1,
      ionization = c("acid", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(4, 0),
      henry_law_constant = 1E-6,
      fu_plasma = 0.5,
      blood_plasma_ratio = 0.8,
      concentration_cells = 2
    )
  )
  expect_error(
    calculate_fu_in_vitro(
      method = "rodgers_rowland_fu",
      lipophilicity = 2,
      ionization = c("base", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22,
      pka = c(6, 0),
      henry_law_constant = 1E-6,
      fu_plasma = 0.3,
      blood_plasma_ratio = 1,
      concentration_cells = 2
    )
  )
})

test_that("calculate_fu_in_vitro: the regressions do not need the well", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "austin",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("base", 0),
      pka = c(3, 0),
      concentration_microsomes = 1
    ),
    0.349395372066127,
    tolerance = 1e-6
  )
  expect_equal(
    calculate_fu_in_vitro(
      method = "kilford",
      system = "cells",
      lipophilicity = 3,
      ionization = c("neutral", "neutral"),
      concentration_cells = 1
    ),
    calculate_fu_in_vitro(
      method = "kilford",
      system = "cells",
      lipophilicity = 3,
      ionization = c("neutral", "neutral"),
      pka = c(0, 0),
      concentration_cells = 1
    )
  )
})

test_that("calculate_fu_in_vitro: berezhkovskiy_fu currently gives the poulin_theil_fu result", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "berezhkovskiy_fu",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "cells",
      fbs_fraction = 0,
      microplate_type = 96,
      fu_plasma = 0.01,
      volume_medium = 0.22,
      pka = c(6, 0),
      henry_law_constant = 1E-6,
      concentration_cells = 2
    ),
    0.783305628878156,
    tolerance = 1e-6
  )
})

test_that("calculate_fu_in_vitro rejects invalid method and system", {
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "bogus",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "cells",
      pka = c(6, 0),
      concentration_cells = 2
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "pksim_standard",
      lipophilicity = 3,
      ionization = c("acid", 0),
      system = "bogus",
      pka = c(6, 0),
      concentration_cells = 2
    )
  )
})

test_that("calculate_fu_in_vitro rejects a method not available for the system", {
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "kilford",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      concentration_microsomes = 1
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "hallifax",
      system = "cells",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      concentration_cells = 1
    )
  )
})

test_that("calculate_fu_in_vitro names the missing arguments", {
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "austin",
      system = "cells",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      concentration_microsomes = 1
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "schmitt_fu",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      concentration_microsomes = 1,
      fbs_fraction = 0,
      microplate_type = 96,
      volume_medium = 0.22
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "pksim_standard",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      concentration_microsomes = 1
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_fu_in_vitro(
      method = "all_literature",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("base", 0),
      pka = c(8, 0),
      concentration_microsomes = 1
    )
  )
})

test_that("calculate_fu_in_vitro warns when the compound probably evaporates", {
  expect_snapshot(
    fu <- calculate_fu_in_vitro(
      method = "austin",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("neutral", 0),
      concentration_microsomes = 1,
      microplate_type = 12,
      volume_medium = 0.5,
      henry_law_constant = 1
    )
  )
})

test_that("calculate_fu_in_vitro: the ionizable group can be in either position", {
  expect_equal(
    calculate_fu_in_vitro(
      method = "austin",
      system = "microsomes",
      lipophilicity = 3,
      ionization = c("neutral", "base"),
      pka = c(0, 8),
      concentration_microsomes = 1
    ),
    0.5689241999032,
    tolerance = 1e-6
  )
})
