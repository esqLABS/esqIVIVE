# Function to convert from Papp to Pint transcellular permeability from a empirical regression

Function to convert from Papp to Pint transcellular permeability from a
empirical regression

## Usage

``` r
pint_peff_empir(Peff_cms, verbose = FALSE)
```

## Arguments

- Peff_cms:

  Peff as obtained in vivo or QSARs in cms

- verbose:

  if TRUE, print the inputs and the resulting pint_cms

## Value

pint in units cm/s

## Examples

``` r
pint_peff_empir(Peff_cms=2.3E-6)
#>     pint_cms 
#> 1.095538e-09 
```
