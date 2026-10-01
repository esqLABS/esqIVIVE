#' Calculate the fraction unbound in an in vitro incubation
#'
#' @description
#' Predicts the fraction unbound of a compound in a microsomal or hepatocyte
#' incubation. Two kinds of method are available:
#'
#' * literature regressions on lipophilicity (`austin`, `hallifax`, `turner`,
#'   `kilford`, `poulin`, and their average `all_literature`), which need only
#'   the compound properties and the microsome or cell concentration;
#' * partition models (`poulin_theil`, `berezhkovskiy`, `pksim_standard`,
#'   `rodgers_rowland`, `schmitt`), which describe binding to the lipids and
#'   proteins of the cells or microsomes, to serum in the medium and to the
#'   plastic of the well. They also need the serum fraction, the microplate
#'   type and the medium volume.
#'
#' @param method Prediction method, one of:
#'
#'   | `method` | System | Description |
#'   |---|---|---|
#'   | `"austin"` | both | Austin et al. (2002) regression |
#'   | `"hallifax"` | microsomes | Hallifax and Houston (2006) regression |
#'   | `"turner"` | microsomes | Turner regression, with separate equations for acids, bases and neutral compounds |
#'   | `"kilford"` | hepatocytes | Kilford et al. (2008) regression |
#'   | `"poulin"` | both | Poulin regression on the neutral lipid content, with acidic phospholipid binding for strong bases |
#'   | `"all_literature"` | both | Average of the regressions available for the system |
#'   | `"poulin_theil"`, `"poulin_theil_fu"` | both | Poulin and Theil partition model |
#'   | `"berezhkovskiy"`, `"berezhkovskiy_fu"` | both | Berezhkovskiy partition model |
#'   | `"pksim_standard"`, `"pksim_standard_fu"` | both | PK-Sim Standard partition model |
#'   | `"rodgers_rowland_fu"` | both | Rodgers and Rowland partition model, strong bases only |
#'   | `"schmitt"`, `"schmitt_fu"` | both | Schmitt partition model |
#'
#'   The partition models without the `_fu` suffix predict binding to serum in
#'   the medium from the serum lipid and protein content. The `_fu` versions
#'   use the measured `fu_plasma` instead.
#' @param system Incubation system, `"microsomes"` or `"hepatocytes"`.
#' @param lipophilicity Lipophilicity of the compound (log units), as logP or
#'   log membrane affinity.
#' @param ionization Ionization class of up to two ionizable groups, as a
#'   vector of length 2 with `"acid"`, `"base"` or `"neutral"`, for example
#'   `c("base", "neutral")`.
#' @param pka pKa values of the two ionizable groups, a vector of length 2.
#'   Defaults to `c(0, 0)`.
#' @param concentration_microsomes Microsomal protein concentration (mg/mL).
#'   Needed when `system = "microsomes"`.
#' @param concentration_cells Hepatocyte concentration (million cells/mL).
#'   Needed when `system = "hepatocytes"`.
#' @param fbs_fraction Fraction of fetal bovine serum in the medium (0 to 1).
#'   Needed by the partition models.
#' @param microplate_type Number of wells of the microplate: 96, 48, 24 or 12.
#'   Needed by the partition models.
#' @param volume_medium Volume of medium in the well (mL). Needed by the
#'   partition models.
#' @param henry_law_constant Henry's law constant (atm m3/mol), used to warn
#'   when the compound probably evaporates from the well. Defaults to 1e-6.
#' @param fu_plasma Fraction unbound in plasma. Needed by the `_fu` methods,
#'   and by `"poulin"` and `"all_literature"` for strong bases.
#' @param blood_plasma_ratio Blood to plasma concentration ratio. Needed by
#'   `"rodgers_rowland_fu"`, and by `"poulin"` and `"all_literature"` for
#'   strong bases.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @return The fraction unbound in the incubation, a single number between 0
#'   and 1. A warning is given when more than 5% of the compound is predicted
#'   to be in the air of the well; this check needs `microplate_type` and
#'   `volume_medium`.
#'
#' @details
#' A strong base is a compound whose first ionization class is `"base"` with
#' a pKa above 7.
#'
#' `"berezhkovskiy"` and `"berezhkovskiy_fu"` currently give the same results
#' as `"poulin_theil"` and `"poulin_theil_fu"`.
#'
#' @concept austin
#' @concept hallifax
#' @concept turner
#' @concept kilford
#' @concept poulin
#' @concept berezhkovskiy
#' @concept rodgers_rowland
#' @concept schmitt
#' @export
#'
#' @examples
#' calculate_fu_in_vitro(
#'   method = "austin",
#'   system = "microsomes",
#'   lipophilicity = 3,
#'   ionization = c("base", "neutral"),
#'   pka = c(8, 0),
#'   concentration_microsomes = 1
#' )
#'
#' calculate_fu_in_vitro(
#'   method = "pksim_standard",
#'   system = "hepatocytes",
#'   lipophilicity = 3,
#'   ionization = c("acid", "neutral"),
#'   pka = c(6, 0),
#'   concentration_cells = 2,
#'   fbs_fraction = 0,
#'   microplate_type = 96,
#'   volume_medium = 0.22
#' )
#'
#' calculate_fu_in_vitro(
#'   method = "poulin_theil_fu",
#'   system = "hepatocytes",
#'   lipophilicity = 3,
#'   ionization = c("acid", "neutral"),
#'   pka = c(6, 0),
#'   concentration_cells = 2,
#'   fbs_fraction = 0.05,
#'   microplate_type = 96,
#'   volume_medium = 0.22,
#'   fu_plasma = 0.01
#' )
calculate_fu_in_vitro <- function(
  method,
  system,
  lipophilicity,
  ionization,
  pka = NULL,
  concentration_microsomes = NULL,
  concentration_cells = NULL,
  fbs_fraction = NULL,
  microplate_type = NULL,
  volume_medium = NULL,
  henry_law_constant = NULL,
  fu_plasma = NULL,
  blood_plasma_ratio = NULL,
  verbose = FALSE
) {
  method <- rlang::arg_match(method, names(.fu_in_vitro_methods))
  system <- rlang::arg_match(system, c("microsomes", "hepatocytes"))
  if (is.null(pka)) {
    pka <- c(0, 0)
  }

  .check_fu_in_vitro_inputs(
    method = method,
    system = system,
    ionization = ionization,
    pka = pka,
    supplied = list(
      concentration_microsomes = concentration_microsomes,
      concentration_cells = concentration_cells,
      fbs_fraction = fbs_fraction,
      microplate_type = microplate_type,
      volume_medium = volume_medium,
      fu_plasma = fu_plasma,
      blood_plasma_ratio = blood_plasma_ratio
    )
  )

  ionization_factors <- .calculate_ionization_factors(ionization, pka)
  ion_factor_plasma <- ionization_factors[["ion_factor_plasma"]] # Interstitial tissue
  ion_factor_cells <- ionization_factors[["ion_factor_cells"]] # intracellular

  protein_partition_1 <- 0.73 * 10^lipophilicity - 0.39 # from Endo 2012,dx.doi.org/10.1021/es303379y partition to chicken muscle, R2=0.86
  protein_partition_2 <- 0.163 + 0.0221 * 10^lipophilicity #from Schmitt 2008, doi:10.1016/j.tiv.2007.09.010
  kPro <- mean(c(protein_partition_1, protein_partition_2))

  # Calculate air-water partition coefficient
  #default to a low hlc if it is not given
  if (is.null(henry_law_constant)) {
    henry_law_constant <- 0.000001
  }
  # Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  # last factor is to convert form atm to Pa
  kAir <- henry_law_constant /
    (0.08206 * 310) *
    101325

  # Calculate plastic partitioning
  plastic_partition_fischer <- 10**(lipophilicity * 0.47 - 4.64)
  plastic_partition_kramer <- 10**(lipophilicity * 0.97 - 6.94)
  kPlastic <- mean(
    c(plastic_partition_fischer, plastic_partition_kramer) *
      1 /
      (1 + ion_factor_cells)
  )

  # get in vitro compartments----------------------------------------------------
  cell <- .calculate_cell_compartments(
    system,
    concentration_cells = concentration_cells,
    concentration_microsomes = concentration_microsomes
  )
  cCellNL <- cell$cell_neutral_lipids
  cCellNPL <- cell$cell_neutral_phospholipids
  cCellAPL <- cell$cell_acidic_phospholipids
  cCellPro <- cell$cell_proteins
  concentration_cell_neutral_lipids <- cCellNL + cCellNPL

  if (.fu_in_vitro_methods[[method]]$partition_model) {
    medium <- .calculate_medium_compartments(fbs_fraction)
    cMediumNL <- medium$medium_neutral_lipids
    cMediumNPL <- medium$medium_neutral_phospholipids
    cMediumPro <- medium$medium_proteins
    well <- .calculate_well_geometry(system, microplate_type, volume_medium)
    saPlasticVolMedium <- well$plastic_area_per_volume
  }

  # QSPRs for calculating partitioning in in vitro------------------------------

  fu_in_vitro <- switch(
    method,
    poulin_theil = ,
    berezhkovskiy = {
      # Calculate lipid partitioning
      kNL <- 10^lipophilicity
      kPL <- 0.3 * 10^lipophilicity + 0.7

      1 /
        (1 +
          kNL * (cCellNL + cMediumNL) +
          kPL * (cCellNPL + cMediumNPL) +
          kPlastic * saPlasticVolMedium)
    },
    poulin_theil_fu = ,
    berezhkovskiy_fu = {
      # Calculate lipid partitioning
      kNL <- 10^lipophilicity
      kNPL <- 0.3 * 10^lipophilicity + 0.7

      1 /
        (1 +
          kNL * cCellNL +
          kNPL * cCellNPL +
          kPlastic * saPlasticVolMedium +
          (1 / fu_plasma - 1) * fbs_fraction)
    },
    pksim_standard = {
      kNL <- 10^lipophilicity

      # assume all neutral lipids have same binding
      1 /
        (1 +
          kNL * (cCellNL + cMediumNL + cCellNPL + cMediumNPL) +
          kPlastic * saPlasticVolMedium +
          kPro * (cCellPro + cMediumPro))
    },
    pksim_standard_fu = {
      kNL <- 10^lipophilicity

      1 /
        (1 +
          kNL * (cCellNL + cCellNPL) +
          kPlastic * saPlasticVolMedium +
          kPro * cCellPro +
          (1 / fu_plasma - 1) * fbs_fraction)
    },
    rodgers_rowland_fu = {
      # RR can only be used by using fu
      # partition into acid phospholipids is only considered if chemical is a strong base
      kOW <- 10^lipophilicity
      kNL <- kOW * (1 / (1 + ion_factor_plasma))
      Hema <- 0.45
      kpuBC <- (Hema - 1 + blood_plasma_ratio) / (Hema * fu_plasma)
      fiwBC <- 0.63
      fnlBC <- 0.003
      fnpBC <- 0.0059
      APbc <- 0.57 # acidic phospholipids in blood cells

      KAPL_1 <- max(
        0,
        kpuBC -
          (1 + ion_factor_cells) / (1 + ion_factor_plasma) * fiwBC -
          (kNL * fnlBC + (0.3 * kNL + 0.7 / (1 + ion_factor_plasma)) * fnpBC)
      )

      kAPL <- KAPL_1 * (1 + ion_factor_plasma) / APbc / ion_factor_cells

      1 /
        (1 +
          kNL * (cCellNL + cMediumNL) +
          (kNL * 0.3 + 0.7 / (1 + ion_factor_plasma)) *
            (cCellNPL + cMediumNPL) +
          kAPL * (cCellAPL) * ion_factor_plasma / (1 + ion_factor_plasma) +
          kPlastic * saPlasticVolMedium)
    },
    schmitt = {
      ionization_parameters_schmitt <- .calculate_ionization_schmitt(
        ionization,
        pka
      )
      logD_Factor <- ionization_parameters_schmitt[["logD_Factor"]]
      kAPLpHFactor <- ionization_parameters_schmitt[["kAPLpHFactor"]]
      LogD <- lipophilicity + log10(logD_Factor)
      kNL <- 10**LogD
      kNPL <- 10**lipophilicity
      kAPL <- kNPL * kAPLpHFactor

      1 /
        (1 +
          kNL * (cCellNL + cMediumNL) +
          kNPL * (cCellNPL + cMediumNPL) +
          kAPL * (cCellAPL) +
          kPro * (cCellPro + cMediumPro))
    },
    schmitt_fu = {
      ionization_parameters_schmitt <- .calculate_ionization_schmitt(
        ionization,
        pka
      )
      logD_Factor <- ionization_parameters_schmitt[["logD_Factor"]]
      kAPLpHFactor <- ionization_parameters_schmitt[["kAPLpHFactor"]]
      LogD <- lipophilicity + log10(logD_Factor)
      kNPL <- 10**lipophilicity
      kNL <- 10**LogD
      kAPL <- kNPL * kAPLpHFactor

      1 /
        (1 +
          kNL * cCellNL +
          kNPL * cCellNPL +
          kAPL * cCellAPL +
          kPro * cCellPro +
          (1 / fu_plasma - 1) * fbs_fraction)
    },
    {
      regressions <- if (method == "all_literature") {
        .fu_in_vitro_methods$all_literature$averages[[system]]
      } else {
        method
      }
      mean(vapply(
        regressions,
        function(regression) {
          .calculate_fu_regression(
            regression,
            system = system,
            ionization = ionization,
            pka = pka,
            lipophilicity = lipophilicity,
            concentration_microsomes = concentration_microsomes,
            concentration_cells = concentration_cells,
            concentration_cell_neutral_lipids = concentration_cell_neutral_lipids,
            concentration_cell_acidic_phospholipids = cCellAPL,
            fu_plasma = fu_plasma,
            blood_plasma_ratio = blood_plasma_ratio
          )
        },
        numeric(1)
      ))
    }
  )

  # Warning for volatility
  if (!is.null(microplate_type) && !is.null(volume_medium)) {
    well <- .calculate_well_geometry(system, microplate_type, volume_medium)
    .check_volatility(fu_in_vitro, kAir, well$volume_air)
  }

  if (verbose) {
    .print_ivive_result(
      "calculate_fu_in_vitro",
      inputs = list(
        method = method,
        system = system,
        lipophilicity = lipophilicity,
        ionization = ionization,
        pka = pka
      ),
      result = fu_in_vitro
    )
  }

  fu_in_vitro
}

# Methods of calculate_fu_in_vitro(): the systems each one is available for,
# whether it is a partition model (which needs the medium and the well), and
# the arguments it needs beyond the concentration for the system.
.fu_in_vitro_partition_arguments <- c(
  "fbs_fraction",
  "microplate_type",
  "volume_medium"
)

.fu_in_vitro_methods <- local({
  regression <- function(systems) {
    list(systems = systems, partition_model = FALSE, requires = character())
  }
  partition <- function(requires = character()) {
    list(
      systems = c("microsomes", "hepatocytes"),
      partition_model = TRUE,
      requires = c(.fu_in_vitro_partition_arguments, requires)
    )
  }
  list(
    austin = regression(c("microsomes", "hepatocytes")),
    hallifax = regression("microsomes"),
    turner = regression("microsomes"),
    kilford = regression("hepatocytes"),
    poulin = regression(c("microsomes", "hepatocytes")),
    all_literature = c(
      regression(c("microsomes", "hepatocytes")),
      list(
        averages = list(
          microsomes = c("poulin", "austin", "hallifax", "turner"),
          hepatocytes = c("kilford", "poulin", "austin")
        )
      )
    ),
    poulin_theil = partition(),
    poulin_theil_fu = partition("fu_plasma"),
    berezhkovskiy = partition(),
    berezhkovskiy_fu = partition("fu_plasma"),
    pksim_standard = partition(),
    pksim_standard_fu = partition("fu_plasma"),
    rodgers_rowland_fu = partition(c("fu_plasma", "blood_plasma_ratio")),
    schmitt = partition(),
    schmitt_fu = partition("fu_plasma")
  )
})

.check_fu_in_vitro_inputs <- function(
  method,
  system,
  ionization,
  pka,
  supplied,
  call = rlang::caller_env()
) {
  definition <- .fu_in_vitro_methods[[method]]
  if (!system %in% definition$systems) {
    cli::cli_abort(
      "{.val {method}} is only available for {definition$systems}.",
      call = call
    )
  }

  concentration <- if (system == "microsomes") {
    "concentration_microsomes"
  } else {
    "concentration_cells"
  }
  requires <- c(concentration, definition$requires)
  missing <- requires[vapply(supplied[requires], is.null, logical(1))]
  if (length(missing) > 0) {
    cli::cli_abort(
      "{.val {method}} needs {.arg {missing}} for {system}.",
      call = call
    )
  }

  # Poulin needs the plasma binding to derive acidic phospholipid binding
  uses_poulin <- method == "poulin" ||
    (method == "all_literature" &&
      "poulin" %in% definition$averages[[system]])
  if (uses_poulin && .is_strong_base(ionization, pka)) {
    requires <- c("fu_plasma", "blood_plasma_ratio")
    missing <- requires[vapply(supplied[requires], is.null, logical(1))]
    if (length(missing) > 0) {
      cli::cli_abort(
        "{.val {method}} needs {.arg {missing}} for strong bases (first
         ionization class {.val base} with a pKa above 7).",
        call = call
      )
    }
  }

  if (method == "rodgers_rowland_fu" && !.is_strong_base(ionization, pka)) {
    cli::cli_abort(
      c(
        "{.val rodgers_rowland_fu} is only available for strong bases: first
       ionization class {.val base} with a pKa above 7.",
        "i" = "PK-Sim uses protein binding instead for acids, neutral compounds
       and weak bases, which is not available here yet."
      ),
      call = call
    )
  }
  invisible(NULL)
}

.is_strong_base <- function(ionization, pka) {
  ionization[1] == "base" && pka[1] > 7
}

.calculate_fu_regression <- function(
  method,
  system,
  ionization,
  pka,
  lipophilicity,
  concentration_microsomes,
  concentration_cells,
  concentration_cell_neutral_lipids,
  concentration_cell_acidic_phospholipids,
  fu_plasma,
  blood_plasma_ratio
) {
  if (method == "poulin") {
    return(.calculate_fu_hep_poulin(
      ionization = ionization,
      pka = pka,
      lipophilicity = lipophilicity,
      concentration_cell_neutral_lipids = concentration_cell_neutral_lipids,
      concentration_cell_acidic_phospholipids = concentration_cell_acidic_phospholipids,
      blood_plasma_ratio = blood_plasma_ratio,
      fu_plasma = fu_plasma
    ))
  }
  if (system == "microsomes") {
    regression <- switch(
      method,
      austin = .calculate_fu_mic_austin,
      hallifax = .calculate_fu_mic_hallifax,
      turner = .calculate_fu_mic_turner
    )
    regression(ionization, pka, lipophilicity, concentration_microsomes)
  } else {
    regression <- switch(
      method,
      austin = .calculate_fu_hep_austin,
      kilford = .calculate_fu_hep_kilford
    )
    regression(ionization, pka, lipophilicity, concentration_cells)
  }
}

.calculate_ionization_schmitt <- function(ionization, pka) {
  # Supports up to 2 ionizable groups (ionization/pka of length 2)
  pH <- 7.4

  ionParam <- c(0, 0)
  for (i in seq(1, 2)) {
    if (ionization[i] == "acid") {
      ionParam[i] <- 1
    } else if (ionization[i] == "base") {
      ionParam[i] <- -1
    } else {
      ionParam[i] <- 0
    }
  }

  # Calculate the fraction neutral
  # conditional if molecule is neutral
  if (abs(ionParam[1]) == 1) {
    F1 <- 1 / (1 + 10^(ionParam[1] * (pka[1] - pH)))
  } else {
    F1 <- 1
  }

  if (abs(ionParam[2]) == 1) {
    F2 <- 1 / (1 + 10^(ionParam[2] * (pka[2] - pH)))
  } else {
    F2 <- 1
  }

  # fraction neutral (both groups uncharged)
  K1 <- F1 * F2
  # fraction with only the first group ionized
  K2 <- (1 - F1) * F2
  # fraction with only the second group ionized
  K3 <- F1 * (1 - F2)
  # fraction with both groups ionized
  K4 <- (1 - F1) * (1 - F2)

  # taken from schmitt paper
  alpha <- 0.001 # ratio of lipophilciity between the neutral and the charged species of a molecule
  # check eq 9 from Schmitt paper
  logD_Factor <- K1 +
    (K2 + K3) * alpha^1 +
    K4 * alpha^max(ionParam[1] + ionParam[2], -ionParam[1] - ionParam[2])

  # check equation 17 and 18 of Schmitt paper
  proportFactorAPL <- 20
  kAPLpHFactor <- K1 +
    K2 * proportFactorAPL^ionParam[1] +
    K3 * proportFactorAPL^ionParam[2] +
    K4 * proportFactorAPL^(ionParam[1] + ionParam[2])

  return(c("logD_Factor" = logD_Factor, "kAPLpHFactor" = kAPLpHFactor))
}


.check_volatility <- function(
  fraction_unbound_in_vitro,
  air_partition_coefficient,
  volume_air_l
) {
  fuAir <- fraction_unbound_in_vitro * air_partition_coefficient * volume_air_l
  # if more than 5 % of the chemicals is predicted to evaporate
  # a warning is given
  # this is conservative because HLC is usually for 25 C and not 37 C
  # and the system is not closed but semi-open,
  # hence a prediction of 5 % with this model actually underpredicts how evaporation will occur

  if (fuAir > 0.05) {
    cli::cli_warn("The compound probably evaporates from the well.")
  }
  invisible(NULL)
}
