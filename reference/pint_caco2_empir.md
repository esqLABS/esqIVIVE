# Function to convert from Papp to Pint transcellular permeability from a empirical regression

Function to convert from Papp to Pint transcellular permeability from a
empirical regression

## Usage

``` r
pint_caco2_empir(Papp_cms, verbose = FALSE)
```

## Arguments

- Papp_cms:

  , permeability from Caco-2 in cms

- verbose:

  if TRUE, print the inputs and the resulting pint_cms

## Value

pint in units cm/s

## Examples

``` r
pint_caco2_empir(Papp_cms=2.3E-6)
#>    pint_cms 
#> 6.16745e-06 
```
