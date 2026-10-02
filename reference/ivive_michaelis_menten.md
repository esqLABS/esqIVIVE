# Scale Michaelis-Menten parameters to in vivo

Scales an in vitro Vmax to the whole tissue with the physiological
scaling factors of the species, and corrects Km for binding in the
incubation.

## Usage

``` r
ivive_michaelis_menten(
  system,
  vmax,
  km,
  fu_in_vitro = 1,
  tissue = "liver",
  species = "human",
  relative_expression_factor = 1,
  verbose = FALSE
)
```

## Arguments

- system:

  Incubation system, `"microsomes"` or `"cells"` (for example
  hepatocytes).

- vmax:

  In vitro Vmax (umol/min/million cells for cells, umol/min/mg protein
  for microsomes), for example from
  [`fit_michaelis_menten_curve()`](https://esqlabs.github.io/ESQivive/reference/fit_michaelis_menten_curve.md).

- km:

  In vitro Km (uM).

- fu_in_vitro:

  Fraction unbound in the incubation, for example from
  [`calculate_fu_in_vitro()`](https://esqlabs.github.io/ESQivive/reference/calculate_fu_in_vitro.md).
  Defaults to 1 (no binding).

- tissue:

  Tissue whose scaling factors are used. Defaults to `"liver"`.

- species:

  Species whose scaling factors are used: `"human"`, `"rat"` or `"dog"`.
  Defaults to `"human"`.

- relative_expression_factor:

  Relative expression or activity factor of the enzyme in vivo compared
  with the incubation. Defaults to 1. To use it, set the reference
  concentration of the enzyme in PK-Sim to 1 uM.

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

A named list:

- `vmax`: in vivo Vmax (umol/min/L of tissue);

- `km_unbound`: unbound Km (uM).

## Examples

``` r
ivive_michaelis_menten(system = "cells", vmax = 2, km = 1)
#> $vmax
#> [1] 355223.9
#> 
#> $km_unbound
#> [1] 1
#> 

ivive_michaelis_menten(
  system = "microsomes",
  vmax = 2,
  km = 1,
  fu_in_vitro = 0.2
)
#> $vmax
#> [1] 107462.7
#> 
#> $km_unbound
#> [1] 0.2
#> 
```
