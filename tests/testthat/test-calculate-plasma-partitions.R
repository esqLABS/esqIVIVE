test_that("predict_plasma_affinities: logP QSAR, ionized compound", {
  # Positional indexing: for an ionized compound, kmemlip_LL/kalb_Lkg inherit a
  # stray name from the ion_factors() lookup (a pre-existing cosmetic quirk,
  # unrelated to the rename fixes made here), so result[["partition_..."]]
  # does not resolve by exact name.
  result <- predict_plasma_affinities(QSAR = "logP", logP = 2, pKa = c(3, 0), ionization = c("acid", 0))

  expect_equal(unname(as.double(result))[1:3], c(1.000183344634427, 0.185104051916421, 0.349000112570181), tolerance = 1e-6)
})

test_that("predict_plasma_affinities: PPLFER QSAR", {
  result <- predict_plasma_affinities(
    QSAR = "PPLFER", logP = 2, pKa = c(3, 0), ionization = c("acid", 0),
    LFER_E = 1, LFER_B = 0, LFER_A = 1.5, LFER_S = 0.8, LFER_V = 2
  )

  expect_equal(
    as.double(result),
    c(11324003.6323556, 16.9112310, 4.8520000),
    tolerance = 1e-6
  )
})

test_that("predict_plasma_affinities rejects an invalid QSAR", {
  expect_error(
    predict_plasma_affinities(QSAR = "bogus", logP = 2, pKa = c(3, 0), ionization = c("acid", 0))
  )
})
