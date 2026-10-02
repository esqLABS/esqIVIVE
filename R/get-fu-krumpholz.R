#' Get measured fraction unbound from the Krumpholz et al. database
#'
#' @description
#' Look up measured fraction unbound values (in microsomes, hepatocytes,
#' plasma or recombinant CYP systems) for one or more compounds in the
#' Krumpholz et al. dataset shipped in `inst/extdata/Krumpholz_et_al_fu_dataset.xlsx`.
#' By default all values measured in the same condition (compound, species
#' and protein or cell concentration) are averaged.
#' Matching of compound names ignores case and leading/trailing spaces. If a
#' compound is not in the database a message is given, with similar compound
#' names when there are any.
#'
#' @param compound character vector with the compound name(s)
#' @param system in vitro system: "microsomes", "hepatocytes", "plasma" or "recombinant CYPs"
#' @param species optional character vector to filter on species (e.g. "human", "rat"). If NULL all species are returned
#' @param average if TRUE (default), average fu over all records of the same
#'  compound, species and concentration. If FALSE, return every record with
#'  method, comments and reference
#' @param verbose if TRUE, print the number of records found per compound
#'
#' @return If `average = TRUE`, a data.frame with one row per condition and columns
#'  compound, species, the concentration of the incubation
#'  (`concentration_mgml` in mg protein/mL for microsomes and recombinant CYPs,
#'  `concentration_Mcellsml` in million cells/mL for hepatocytes, none for plasma),
#'  `cyp` (recombinant CYPs only), `fu` (mean) and `n` (number of values averaged).
#'  Records whose concentration was reported as a range or not reported are
#'  averaged together under an NA concentration. Duplicate records (same
#'  compound, species, concentration, fu and reference) are counted once.
#'
#'  If `average = FALSE`, a data.frame with one row per record and columns
#'  compound, system, species, method, concentration, concentration_unit,
#'  concentration_reported, fu, fu_sd, fu_range, comments, doi and reference.
#'  For "recombinant CYPs" the columns cyp, concentration_pmol_cyp_ml and
#'  compound_concentration_uM are added.
#'
#'  Compounds that are not found return no rows.
#' @export
#'
#' @examples
#' # all verapamil fu_mic in human, averaged per microsomal concentration
#' get_fu_krumpholz("Verapamil", system = "microsomes", species = "human")
#'
#' get_fu_krumpholz(c("midazolam", "Diazepam"), system = "hepatocytes")
#'
#' # individual records with references
#' get_fu_krumpholz("Verapamil", system = "microsomes", species = "human", average = FALSE)
#'
#' # check if a compound is in the database
#' nrow(get_fu_krumpholz("Verapamil", system = "plasma")) > 0
get_fu_krumpholz <- function(
  compound,
  system = "microsomes",
  species = NULL,
  average = TRUE,
  verbose = FALSE
) {
  rlang::arg_match(system, names(.krumpholz_sheets))

  fu_data <- .read_krumpholz(system)
  compound_key <- .normalize_name(compound)
  data_key <- .normalize_name(fu_data$compound)

  missing <- compound[!compound_key %in% data_key]
  for (cmp in missing) {
    similar <- unique(fu_data$compound[agrepl(cmp, fu_data$compound, max.distance = 0.2, ignore.case = TRUE)])
    message(
      "\"", cmp, "\" not found in the Krumpholz ", system, " data.",
      if (length(similar) > 0) {
        paste0(" Did you mean: ", paste(utils::head(similar, 5), collapse = ", "), "?")
      }
    )
  }

  result <- fu_data[data_key %in% compound_key, ]
  if (!is.null(species)) {
    result <- result[tolower(result$species) %in% tolower(species), ]
  }
  rownames(result) <- NULL

  if (verbose) {
    found <- table(factor(.normalize_name(result$compound), levels = unique(compound_key)))
    .print_ivive_result(
      "get_fu_krumpholz",
      inputs = list(compound = compound, system = system, species = species),
      result = stats::setNames(as.integer(found), names(found))
    )
  }

  if (average) {
    result <- .average_fu(result, system)
  }

  result
}

# Mean fu per compound, species and concentration (and CYP for recombinant CYPs)
.average_fu <- function(fu_data, system) {
  fu_data <- fu_data[!is.na(fu_data$fu), ]
  group_columns <- c("compound", "species", "concentration")
  if (system == "recombinant CYPs") {
    group_columns <- c(group_columns, "cyp")
  }

  # the same value from the same reference in the same condition is only
  # counted once (the database has repeated entries, e.g. with and without comments)
  fu_data <- fu_data[!duplicated(fu_data[c(group_columns, "concentration_reported", "fu", "reference")]), ]

  if (nrow(fu_data) == 0) {
    averaged <- data.frame(fu_data[group_columns], fu = numeric(0), n = integer(0))
  } else {
    # NA is kept as its own group (e.g. concentration reported as a range)
    group <- interaction(lapply(fu_data[group_columns], addNA), drop = TRUE, lex.order = TRUE)
    averaged <- do.call(rbind, lapply(split(fu_data, group), function(x) {
      data.frame(x[1, group_columns, drop = FALSE], fu = mean(x$fu), n = nrow(x))
    }))
    averaged <- averaged[order(averaged$compound, averaged$species, averaged$concentration), ]
  }
  rownames(averaged) <- NULL

  concentration_name <- switch(
    system,
    "hepatocytes" = "concentration_Mcellsml",
    "plasma" = NA,
    "concentration_mgml"
  )
  if (is.na(concentration_name)) {
    averaged$concentration <- NULL
  } else {
    names(averaged)[names(averaged) == "concentration"] <- concentration_name
  }

  averaged
}

#' List compounds in the Krumpholz et al. database
#'
#' @param system in vitro system: "microsomes", "hepatocytes", "plasma" or "recombinant CYPs"
#'
#' @return sorted character vector with the compound names available for that system
#' @export
#'
#' @examples
#' head(list_krumpholz_compounds("hepatocytes"))
list_krumpholz_compounds <- function(system = "microsomes") {
  rlang::arg_match(system, names(.krumpholz_sheets))
  sort(unique(.read_krumpholz(system)$compound))
}

# sheet name and concentration column/unit per system
.krumpholz_sheets <- list(
  "microsomes" = list(sheet = "Microsomes", conc = "mg/ml", unit = "mg protein/mL"),
  "hepatocytes" = list(sheet = "Hepatocytes", conc = "x10^6 cells", unit = "million cells/mL"),
  "plasma" = list(sheet = "Plasma", conc = NA, unit = NA),
  "recombinant CYPs" = list(sheet = "Recombinant CYPs", conc = "concentration [mg protein/mL]", unit = "mg protein/mL")
)

# cache of the sheets already read in this session
.krumpholz_cache <- new.env(parent = emptyenv())

.read_krumpholz <- function(system) {
  if (!is.null(.krumpholz_cache[[system]])) {
    return(.krumpholz_cache[[system]])
  }

  info <- .krumpholz_sheets[[system]]
  path <- system.file("extdata", "Krumpholz_et_al_fu_dataset.xlsx", package = "ESQivive")
  raw <- as.data.frame(readxl::read_excel(path, sheet = info$sheet, col_types = "text"))

  column <- function(name) {
    if (is.na(name) || !name %in% names(raw)) rep(NA_character_, nrow(raw)) else raw[[name]]
  }
  # "0.0008*" -> 0.0008; ranges ("0.5-1") and codes ("ND", "NC") -> NA
  as_number <- function(x) suppressWarnings(as.numeric(sub("*", "", x, fixed = TRUE)))

  concentration_reported <- column(info$conc)
  fu_data <- data.frame(
    compound = trimws(raw$Compound),
    system = system,
    species = column("species"),
    method = column("method"),
    concentration = as_number(concentration_reported),
    concentration_unit = info$unit,
    concentration_reported = concentration_reported,
    fu = as_number(raw[["fu mean"]]),
    fu_sd = as_number(raw$SD),
    fu_range = column("fu range"),
    comments = column("comments"),
    doi = column("REF"),
    reference = column("Author, year")
  )

  if (system == "recombinant CYPs") {
    fu_data$cyp <- column("CYP")
    fu_data$concentration_pmol_cyp_ml <- as_number(column("concentration [pmol CYP/mL]"))
    fu_data$compound_concentration_uM <- as_number(column("concentration [uM]"))
  }

  .krumpholz_cache[[system]] <- fu_data
  fu_data
}

.normalize_name <- function(x) tolower(trimws(x))
