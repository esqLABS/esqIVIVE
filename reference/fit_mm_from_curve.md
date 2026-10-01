# Get michaelis-menten parameters form experimental curves

Function to derive micahelis-menten form raw data

## Usage

``` r
fit_mm_from_curve(experimental_conc_velocity, verbose = FALSE)
```

## Arguments

- experimental_conc_velocity:

  is a experimental curve with concentration in the first column and
  velocity in the second column

- verbose:

  if TRUE, print the inputs and the resulting Km/Vmax

## Value

fitresults_vmax_km

## Examples

``` r
mm_curve_path<-system.file("extdata","michaelis_menten_curve.csv",package="esqIVIVE")
mm_curve<-read.csv(mm_curve_path)
fit_mm_from_curve(mm_curve)
#> Waiting for profiling to be done...

#>                                    Mean X2.5_percent X95._percent
#> Km_uM                        24.8194106   17.5269630   35.5214188
#> Vmax_umol_min_mgmicroORcells  0.1480423    0.1307565    0.1700504
```
