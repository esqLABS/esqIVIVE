# ion_factors

Helper Function to calculate fraction unionized If there are multiple
pKas for acidity just use the lower value If there are multiple pKas for
basicity just use the lower value Mind that pKb is not the same as pKa !
If you have pKb, just calculate pKa=14-pKb

## Usage

``` r
ion_factors(ionization, pKa, verbose = FALSE)
```

## Arguments

- ionization:

  vector of length 2 with ionization class, acid, neutral and base

- pKa:

  vector of length 2 with pKa values of the compound

- verbose:

  if TRUE, print the inputs and the resulting ionization factors

## Value

factors that can be used to calculate the fraction neutral or ionized in
plasma and intracellularly

## Examples

``` r
ion_factors(ionization=c("neutral",0),pKa<-c(0,0))
#> ion_factor_plasma  ion_factor_cells 
#>                 0                 0 
ion_factors(ionization=c("acid",0),pKa<-c(14,0))
#> ion_factor_plasma  ion_factor_cells 
#>      2.511886e-07      1.659587e-07 
ion_factors(ionization=c("base","acid"),pKa<-c(5,7))
#> ion_factor_plasma  ion_factor_cells 
#>          2.515868          1.665613 
```
