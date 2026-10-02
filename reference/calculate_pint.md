# Calculate the intestinal transcellular permeability

Converts a measured permeability into the intestinal transcellular
permeability (Pint) used by PK-Sim, with an empirical regression
calibrated on compounds of high solubility.

## Usage

``` r
calculate_pint(method, permeability, verbose = FALSE)
```

## Arguments

- method:

  Type of the measured permeability: `"caco2"` for the apparent
  permeability in Caco-2 cells (Papp), or `"peff"` for the human
  effective intestinal permeability (Peff), measured in vivo or
  predicted.

- permeability:

  Measured permeability (cm/s).

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

The intestinal transcellular permeability (cm/s), a single number.

## Examples

``` r
calculate_pint(method = "caco2", permeability = 2.3e-6)
#> [1] 6.16745e-06

calculate_pint(method = "peff", permeability = 2.3e-6)
#> [1] 1.095538e-09
```
