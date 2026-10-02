# Calculate the compartments of an in vitro incubation

Describes an incubation of cells (such as hepatocytes) or microsomes in
a microplate well: the lipid and protein content of the cells or
microsomes and of the serum in the medium, the plastic surface in
contact with the medium, and the air above it.
[`calculate_fu_in_vitro()`](https://esqlabs.github.io/ESQivive/reference/calculate_fu_in_vitro.md)
uses these values for its partition models. You can also use them to
describe a virtual incubation.

## Usage

``` r
calculate_in_vitro_compartments(
  system,
  fbs_fraction,
  microplate_type,
  volume_medium,
  concentration_cells = NULL,
  concentration_microsomes = NULL,
  verbose = FALSE
)
```

## Arguments

- system:

  Incubation system, `"microsomes"` or `"cells"` (for example
  hepatocytes).

- fbs_fraction:

  Fraction of fetal bovine serum in the medium (0 to 1).

- microplate_type:

  Number of wells of the microplate: 96, 48, 24 or 12.

- volume_medium:

  Volume of medium in the well (mL).

- concentration_cells:

  Cell concentration (million cells/mL). Needed when `system = "cells"`.

- concentration_microsomes:

  Microsomal protein concentration (mg/mL). Needed when
  `system = "microsomes"`.

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

A named list. Lipid and protein contents are fractions of the medium
volume (L/L).

- `cell_neutral_lipids`, `cell_neutral_phospholipids`,
  `cell_acidic_phospholipids`, `cell_proteins`: content of the cells or
  microsomes.

- `medium_neutral_lipids`, `medium_neutral_phospholipids`,
  `medium_proteins`: content of the serum in the medium.

- `plastic_area_per_volume`: plastic surface in contact with the medium
  per medium volume (m2/L). Zero for microsomes, which are assumed to be
  incubated in glass.

- `volume_air`: volume of air above the medium in the well (L).

## Examples

``` r
calculate_in_vitro_compartments(
  system = "cells",
  fbs_fraction = 0.05,
  microplate_type = 96,
  volume_medium = 0.15,
  concentration_cells = 0.1
)
#> $cell_neutral_lipids
#> [1] 1.1303e-05
#> 
#> $cell_neutral_phospholipids
#> [1] 8.4074e-06
#> 
#> $cell_acidic_phospholipids
#> [1] 2.2352e-06
#> 
#> $cell_proteins
#> [1] 5.08e-05
#> 
#> $medium_neutral_lipids
#> [1] 7.85e-05
#> 
#> $medium_neutral_phospholipids
#> [1] 1.5e-05
#> 
#> $medium_proteins
#> [1] 0.002
#> 
#> $plastic_area_per_volume
#> [1] 0.6060606
#> 
#> $volume_air
#> [1] 0.000242
#> 

calculate_in_vitro_compartments(
  system = "microsomes",
  fbs_fraction = 0,
  microplate_type = 24,
  volume_medium = 0.5,
  concentration_microsomes = 1
)
#> $cell_neutral_lipids
#> [1] 0.0002611111
#> 
#> $cell_neutral_phospholipids
#> [1] 0.0007261556
#> 
#> $cell_acidic_phospholipids
#> [1] 0.0001594
#> 
#> $cell_proteins
#> [1] 0.0007407407
#> 
#> $medium_neutral_lipids
#> [1] 0
#> 
#> $medium_neutral_phospholipids
#> [1] 0
#> 
#> $medium_proteins
#> [1] 0
#> 
#> $plastic_area_per_volume
#> [1] 0
#> 
#> $volume_air
#> [1] 0.00297
#> 
```
