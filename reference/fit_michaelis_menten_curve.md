# Fit a Michaelis-Menten curve

Fits the Michaelis-Menten equation to reaction velocities measured at
several substrate concentrations and returns Km and Vmax. A plot of the
data and the fitted curve is drawn so you can judge the fit.

To scale the results to in vivo, pass them to
[`ivive_michaelis_menten()`](https://esqlabs.github.io/ESQivive/reference/ivive_michaelis_menten.md).

## Usage

``` r
fit_michaelis_menten_curve(data, verbose = FALSE)
```

## Arguments

- data:

  A data frame with the substrate concentration (uM) in the first column
  and the velocity in the second column. Rows with a missing
  concentration or velocity are left out.

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

A data frame with two rows, `parameter = "km"` (uM) and
`parameter = "vmax"` (in the velocity unit of `data`), and the columns
`estimate`, `lower` and `upper`, the bounds of the 95% confidence
interval. A warning is given when the fit is poor (R-squared below 0.8).

## Examples

``` r
mm_curve <- read.csv(
  system.file("extdata", "michaelis_menten_curve.csv", package = "ESQivive")
)
fit_michaelis_menten_curve(mm_curve)
#> Waiting for profiling to be done...

#>   parameter   estimate      lower      upper
#> 1        km 24.8194106 17.5269630 35.5214188
#> 2      vmax  0.1480423  0.1307565  0.1700504
```
