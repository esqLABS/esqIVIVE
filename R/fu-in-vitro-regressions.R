# Literature regressions for the fraction unbound in microsomal and hepatocyte
# incubations. Internal: reached through calculate_fu_in_vitro(), which checks
# the inputs. Concentrations are in mg protein/mL (microsomes) and million
# cells/mL (hepatocytes).

# Lipophilicity used by the regressions: logD for strong bases, logP otherwise.
.regression_lipophilicity <- function(ionization, pka, lipophilicity) {
  if (.is_strong_base(ionization, pka)) {
    ion_factor <- .calculate_ionization_factors(ionization, pka)[[
      "ion_factor_plasma"
    ]]
    log10(1 / (1 + ion_factor) * 10^lipophilicity)
  } else {
    lipophilicity
  }
}

# Turner et al 2006 (Turner-Simcyp model), Drug Metab Rev 38(S1):162, a
# conference abstract without a DOI. See the references in calculate_fu_in_vitro().
.calculate_fu_mic_turner <- function(
  ionization,
  pka,
  lipophilicity,
  concentration_microsomes
) {
  if (.is_strong_base(ionization, pka)) {
    1 / (1 + concentration_microsomes * 10^(0.58 * lipophilicity - 2.02))
  } else if (ionization[1] == "acid" && pka[1] < 7) {
    1 / (1 + concentration_microsomes * 10^(0.2 * lipophilicity - 1.54))
  } else {
    1 / (1 + concentration_microsomes * 10^(0.46 * lipophilicity - 1.51))
  }
}

# Hallifax and Houston 2006, doi:10.1124/dmd.105.007658
.calculate_fu_mic_hallifax <- function(
  ionization,
  pka,
  lipophilicity,
  concentration_microsomes
) {
  log_partition <- .regression_lipophilicity(ionization, pka, lipophilicity)
  1 /
    (1 +
      concentration_microsomes *
        10^(0.072 * log_partition^2 + 0.067 * log_partition - 1.126))
}

# Austin et al 2002, doi:10.1124/dmd.30.12.1497
.calculate_fu_mic_austin <- function(
  ionization,
  pka,
  lipophilicity,
  concentration_microsomes
) {
  log_partition <- .regression_lipophilicity(ionization, pka, lipophilicity)
  1 / (1 + concentration_microsomes * 10^(0.56 * log_partition - 1.41))
}

# Austin et al 2005, doi:10.1124/dmd.104.002436: the hepatocyte regression
# log((1 - fu)/fu) = 0.4 log(P or D) - 1.38 per million cells/mL
.calculate_fu_hep_austin <- function(
  ionization,
  pka,
  lipophilicity,
  concentration_cells
) {
  log_partition <- .regression_lipophilicity(ionization, pka, lipophilicity)
  1 / (1 + concentration_cells * 10^(0.4 * log_partition - 1.38))
}

# Kilford et al 2008, doi:10.1124/dmd.108.020834
.calculate_fu_hep_kilford <- function(
  ionization,
  pka,
  lipophilicity,
  concentration_cells
) {
  log_partition <- .regression_lipophilicity(ionization, pka, lipophilicity)
  volume_ratio <- 0.005 * concentration_cells
  1 /
    (1 +
      125 *
        volume_ratio *
        10^(0.072 * log_partition^2 + 0.067 * log_partition - 1.126))
}

# Poulin: binding to the neutral lipids of the cells or microsomes, plus the
# acidic phospholipids for strong bases. Used for both systems. Lipid
# concentrations are fractions of the medium volume.
# Poulin and Haddad 2011 (microsomes), doi:10.1002/jps.22619, and 2013 (cells),
# doi:10.1002/jps.23602. The acidic phospholipid partition is calibrated on
# blood cells, following Rodgers et al 2005, doi:10.1002/jps.20322.
.calculate_fu_hep_poulin <- function(
  ionization,
  pka,
  lipophilicity,
  concentration_cell_neutral_lipids,
  concentration_cell_acidic_phospholipids = NULL,
  blood_plasma_ratio = NULL,
  fu_plasma = NULL
) {
  neutral_lipid_partition <- 10^lipophilicity
  ionization_factors <- .calculate_ionization_factors(ionization, pka)
  ion_factor_plasma <- ionization_factors[["ion_factor_plasma"]]
  # the acidic phospholipid partition is calibrated on blood cells (Rodgers and Rowland)
  ion_factor_blood_cells <- ionization_factors[["ion_factor_blood_cells"]]

  if (.is_strong_base(ionization, pka)) {
    fraction_neutral_lipids_erythrocytes <- 0.0024
    fraction_acidic_phospholipids_erythrocytes <- 0.00057
    fraction_water_erythrocytes <- 0.63
    partition_erythrocytes_albumin <- (blood_plasma_ratio - (1 - 0.45)) /
      0.45 /
      fu_plasma

    acidic_phospholipid_partition <- (partition_erythrocytes_albumin -
      ((1 + ion_factor_blood_cells) *
        fraction_water_erythrocytes +
        neutral_lipid_partition * fraction_neutral_lipids_erythrocytes) /
        (1 + ion_factor_plasma)) *
      ((1 + ion_factor_plasma) /
        (ion_factor_blood_cells * fraction_acidic_phospholipids_erythrocytes))

    1 /
      (1 +
        ((neutral_lipid_partition *
          concentration_cell_neutral_lipids +
          ion_factor_plasma *
            acidic_phospholipid_partition *
            concentration_cell_acidic_phospholipids) /
          (1 + ion_factor_plasma)))
  } else {
    1 /
      (1 +
        ((neutral_lipid_partition * concentration_cell_neutral_lipids) /
          (1 + ion_factor_plasma)))
  }
}
