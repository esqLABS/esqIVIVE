# Correct the fraction unbound in plasma for binding to neutral lipids

Applies the Pearce correction to a measured fraction unbound in plasma,
to account for binding to the neutral lipids of plasma that is not seen
in the measurement.

## Usage

``` r
correct_fu_plasma_pearce(fu_plasma, lipophilicity, verbose = FALSE)
```

## Arguments

- fu_plasma:

  Measured fraction unbound in plasma.

- lipophilicity:

  Lipophilicity of the compound (log units), as logP or log membrane
  affinity.

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

The corrected fraction unbound in plasma, a single number.

## Examples

``` r
correct_fu_plasma_pearce(fu_plasma = 0.2, lipophilicity = 4)
#> [1] 0.01333333
```
