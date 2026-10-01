# Correct Fu based on Pearce correction factor

Corrects the Fu based on Pearce correction factor for neutral lipids in
plasma

## Usage

``` r
correct_fu_pls_pearce(fraction_unbound, log_lipophilicity, verbose = FALSE)
```

## Arguments

- fraction_unbound:

  Fraction unbound in plasma

- log_lipophilicity:

  LogP or LogMA of the compound

- verbose:

  if TRUE, print the inputs and the resulting Fu_plasma

## Value

Corrected Fu_plasma value

## Examples

``` r
correct_fu_pls_pearce(fraction_unbound=0.2, log_lipophilicity=4)
#>  Fu_plasma 
#> 0.01333333 
```
