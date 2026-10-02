test_that("calculate_plasma_partitions: logp method, ionized compound", {
  # Tests the logP-based method for an ionized compound. A strong acid
  # (pKa 3) is almost fully ionized at plasma pH 7.4, so its lipophilicity is
  # first corrected to logD and the partition coefficients come out very low.
  # The expected values follow from the equations in the method:
  # kmemlip = 10^logD = 10^2 / (1 + 10^(7.4 - 3)); kalb = 0.163 + 0.0221 * kmemlip;
  # kglob = mean(kalb, 10^(0.37 * logD - 0.29))
  expect_equal(
    calculate_plasma_partitions(
      method = "logp",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0)
    ),
    list(
      partition_albumin = 0.16308797818221782,
      partition_globulin = 0.11473065377891481,
      partition_membrane_lipids = 0.00398091322252505
    ),
    tolerance = 1e-6
  )
})

test_that("calculate_plasma_partitions: pplfer method matches the LSER Calculation workbook", {
  # Tests the PP-LFER method against reference values from the "LSER
  # Calculation" workbook, for four chemicals (ketoconazole, diclofenac,
  # naproxen and diphenhydramine). The workbook is independent of the
  # coefficients in the code, so it checks the equations were entered correctly.
  #
  # For each chemical the workbook gives the Abraham descriptors and the logK
  # of each system from equation 1 (L-based) and equation 3 (E-based), rounded
  # to 2 decimals. Expected K = mean of the two K values (linear space); the
  # protein systems are converted from L/L to L/kg with a protein density of
  # 1.35 g/mL.
  chemicals <- list(
    ketoconazole = list(
      d = c(e = 3.14, s = 3.68, a = 0, b = 2.6, v = 3.7208, l = 19.71),
      muscle = c(3.08, 2.53),
      membrane = c(3.45, 2.80),
      albumin = c(2.88, 2.54)
    ),
    diclofenac = list(
      d = c(e = 1.81, s = 1.85, a = 0.55, b = 0.77, v = 2.025, l = 11.025),
      muscle = c(3.59, 3.27),
      membrane = c(4.73, 4.25),
      albumin = c(4.10, 3.87)
    ),
    naproxen = list(
      d = c(e = 1.51, s = 2.022, a = 0.6, b = 0.673, v = 1.7821, l = 9.207),
      muscle = c(2.69, 2.60),
      membrane = c(3.61, 3.46),
      albumin = c(3.39, 3.36)
    ),
    diphenhydramine = list(
      d = c(e = 1.31, s = 1.11, a = 0, b = 1.22, v = 2.1872, l = 9.31),
      muscle = c(2.27, 2.40),
      membrane = c(3.27, 3.25),
      albumin = c(2.72, 2.68)
    )
  )

  for (chemical in names(chemicals)) {
    x <- chemicals[[chemical]]
    result <- calculate_plasma_partitions(
      method = "pplfer",
      lfer_e = x$d[["e"]],
      lfer_s = x$d[["s"]],
      lfer_a = x$d[["a"]],
      lfer_b = x$d[["b"]],
      lfer_v = x$d[["v"]],
      lfer_l = x$d[["l"]]
    )

    expect_named(
      result,
      c("partition_albumin", "partition_globulin", "partition_membrane_lipids")
    )
    # the workbook rounds logK to 2 decimals, i.e. up to ~1.2% in K
    expect_equal(
      result$partition_membrane_lipids,
      mean(10^x$membrane),
      tolerance = 0.015,
      info = chemical
    )
    expect_equal(
      result$partition_albumin,
      mean(10^x$albumin) / 1.35,
      tolerance = 0.015,
      info = chemical
    )
    expect_equal(
      result$partition_globulin,
      mean(10^x$muscle) / 1.35,
      tolerance = 0.015,
      info = chemical
    )
  }
})

test_that("calculate_plasma_partitions: pplfer method ignores ionization inputs", {
  # Tests that the PP-LFER equations, which are for the neutral species, give
  # the same result whether or not lipophilicity, ionization and pKa are
  # supplied, and that the method can be called without them at all.
  descriptors <- list(
    lfer_e = 0.610,
    lfer_s = 0.52,
    lfer_a = 0,
    lfer_b = 0.14,
    lfer_v = 0.7164,
    lfer_l = 2.786
  )

  expect_equal(
    do.call(
      calculate_plasma_partitions,
      c(
        list(
          method = "pplfer",
          lipophilicity = 2,
          pka = c(3, 0),
          ionization = c("acid", 0)
        ),
        descriptors
      )
    ),
    do.call(calculate_plasma_partitions, c(list(method = "pplfer"), descriptors))
  )
})

test_that("calculate_plasma_partitions rejects invalid inputs", {
  # Tests the error messages, stored as snapshots, for an unknown method and
  # for the PP-LFER method called with missing descriptors (the message must
  # list every missing descriptor, including lfer_l).
  expect_snapshot(
    error = TRUE,
    calculate_plasma_partitions(
      method = "bogus",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0)
    )
  )
  expect_snapshot(
    error = TRUE,
    calculate_plasma_partitions(
      method = "pplfer",
      lipophilicity = 2,
      pka = c(3, 0),
      ionization = c("acid", 0),
      lfer_e = 1
    )
  )
})

test_that("calculate_plasma_partitions: the ionizable group can be in either position", {
  # Tests that a single ionizable group gives the same result in the second
  # slot (c("neutral", "acid"), pKa c(0, 3)) as in the first (c("acid",
  # "neutral"), pKa c(3, 0)), which is the case in the test above.
  expect_equal(
    calculate_plasma_partitions(
      method = "logp",
      lipophilicity = 2,
      pka = c(0, 3),
      ionization = c("neutral", "acid")
    ),
    list(
      partition_albumin = 0.16308797818221782,
      partition_globulin = 0.11473065377891481,
      partition_membrane_lipids = 0.00398091322252505
    ),
    tolerance = 1e-6
  )
})
