# Calculate fraction unbound in the in vitro (hepatic) system

Compute the fraction unbound in vitro with the option of using different
QSARs

## Usage

``` r
calculate_fu_in_vitro(
  partition_qspr,
  log_lipophilicity,
  ionization,
  type_system,
  FBS_fraction,
  microplate_type,
  volume_medium,
  pka = NULL,
  henry_law_constant = NULL,
  fraction_unbound = NULL,
  blood_plasma_ratio = NULL,
  concentration_microsomes = NULL,
  concentration_cells = NULL,
  verbose = FALSE
)
```

## Arguments

- partition_qspr:

  type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers &
  Rowland, Schmidtt, then from literature, Poulin, Turner, Austin and
  Halifax. See thevignette for more details)

- log_lipophilicity:

  LogP or LogMA of the compound

- ionization:

  Vector of length 2 with ionization class, acid, neutral and base, if
  not input then it is c(0,0)

- type_system:

  microsomes or hepatocytes

- FBS_fraction:

  fraction of serum concentration, values can only go from 0-1

- microplate_type:

  number of wells in the microplate

- volume_medium:

  volume of medium in the well (in mL)

- pka:

  vector of length of 2 with pkA of the compound

- henry_law_constant:

  Henry's Law Constant (in atm/(m3\*mol))

- fraction_unbound:

  In Vivo Fraction Unbound in plasma from literature

- blood_plasma_ratio:

  Blood plasma ratio, this parameter is needed for Rodgers and Rowland
  and Poulin method for basic chemicals

- concentration_microsomes:

  concentration of microsomes (in mg/mL)

- concentration_cells:

  concentration of cells (in million cells/mL)

- verbose:

  if TRUE, print the inputs and the resulting fuInvitro

## Value

fuInvitro and possible warning for evaporation

## Details

mayeb consider to have average data..

## Examples

``` r
calculate_fu_in_vitro(
 partition_qspr = "All PK-Sim Standard", log_lipophilicity = 3, ionization = c("acid", 0),
 type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
 volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2)
#> [1] 0.5630088

calculate_fu_in_vitro(
 partition_qspr = "Poulin and Theil + fu", log_lipophilicity = 3, ionization = c("acid", 0),
 type_system = "hepatocytes", FBS_fraction = 0, microplate_type = 96,
 fraction_unbound=0.01,blood_plasma_ratio=2,
 volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2)
#> [1] 0.78331

calculate_fu_in_vitro(
 partition_qspr = "All Schmitt", log_lipophilicity = 0.42, ionization = c("acid", 0),
 type_system = "microsomes", FBS_fraction = 0, microplate_type = 96,
 fraction_unbound=0.2,blood_plasma_ratio=1,
 volume_medium = 0.22, pka = c(6, 0), concentration_microsomes = 2)
#> [1] 0.992175
```
