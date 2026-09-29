test_that("get_fu_krumpholz: averages fu per compound, species and concentration", {
  result <- get_fu_krumpholz("  verapamil ", system = "microsomes", species = "human")
  expect_named(result, c("compound", "species", "concentration_mgml", "fu", "n"))
  expect_true(all(result$compound == "Verapamil"))
  expect_true(all(result$species == "human"))
  # one row per concentration (NA included once)
  expect_equal(anyDuplicated(result$concentration_mgml), 0)

  records <- get_fu_krumpholz("Verapamil", system = "microsomes", species = "human", average = FALSE)
  at_1 <- unique(records[records$concentration %in% 1, c("fu", "reference")])
  expect_equal(result$fu[result$concentration_mgml %in% 1], mean(at_1$fu))
  expect_equal(result$n[result$concentration_mgml %in% 1], nrow(at_1))
})

test_that("get_fu_krumpholz: duplicate records are counted once in the average", {
  records <- get_fu_krumpholz("Verapamil", system = "microsomes", species = "human", average = FALSE)
  mclure <- records[records$reference == "McLure 2011" & records$fu == 0.7, ]
  expect_equal(nrow(mclure), 2)

  result <- get_fu_krumpholz("Verapamil", system = "microsomes", species = "human")
  at_1 <- records[records$concentration %in% 1, ]
  expect_equal(result$n[result$concentration_mgml %in% 1], nrow(at_1) - 1)
})

test_that("get_fu_krumpholz: concentrations reported as a range are grouped as NA", {
  result <- get_fu_krumpholz("Verapamil", system = "microsomes", species = "human")
  expect_equal(sum(is.na(result$concentration_mgml)), 1)
})

test_that("get_fu_krumpholz: hepatocytes, plasma and several compounds", {
  hep <- get_fu_krumpholz(c("Acetaminophen", "Albendazole"), system = "hepatocytes")
  expect_setequal(unique(hep$compound), c("Acetaminophen", "Albendazole"))
  expect_true("concentration_Mcellsml" %in% names(hep))

  pls <- get_fu_krumpholz("Acetaminophen", system = "plasma")
  expect_named(pls, c("compound", "species", "fu", "n"))
})

test_that("get_fu_krumpholz: individual records with average = FALSE", {
  result <- get_fu_krumpholz("Verapamil", system = "microsomes", species = "human", average = FALSE)
  expect_true(all(result$concentration_unit == "mg protein/mL"))
  obach <- result[result$reference == "Obach 1999", ]
  expect_equal(obach$fu, 0.43)
  expect_equal(obach$concentration, 0.5)
  range_row <- result[which(result$concentration_reported == "0.5-1"), ]
  expect_true(is.na(range_row$concentration))
})

test_that("get_fu_krumpholz: unknown compound gives a message and no rows", {
  expect_message(
    result <- get_fu_krumpholz("Verapamill", system = "microsomes"),
    "Did you mean: Verapamil"
  )
  expect_equal(nrow(result), 0)
  expect_named(result, c("compound", "species", "concentration_mgml", "fu", "n"))
})

test_that("get_fu_krumpholz: recombinant CYPs are also grouped by CYP", {
  result <- get_fu_krumpholz("Alprazolam", system = "recombinant CYPs")
  expect_true(all(c("cyp", "concentration_mgml") %in% names(result)))
  expect_equal(result$cyp, "3A4")
})

test_that("get_fu_krumpholz: invalid system errors", {
  expect_error(get_fu_krumpholz("Verapamil", system = "kidney"))
})

test_that("list_krumpholz_compounds: lists available compounds", {
  compounds <- list_krumpholz_compounds("microsomes")
  expect_true("Verapamil" %in% compounds)
  expect_false(is.unsorted(compounds))
})
