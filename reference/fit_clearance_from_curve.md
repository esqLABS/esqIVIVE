# function to determine clearance from experimental curve

function to determine clearance from experimental curve

## Usage

``` r
fit_clearance_from_curve(expData_tmin_cuM, verbose = FALSE)
```

## Arguments

- expData_tmin_cuM:

  table with first column as time in minutes and second column as
  concentration in uM

- verbose:

  if TRUE, print the inputs and the resulting kcat

## Value

kcat in per min, this is not clearance ready for pksim

## Examples

``` r
exp_path<-system.file("extdata","clearance.csv",package="esqIVIVE")
expData<-read.csv(exp_path)
#see that the first column is time and second is concentration
expData
#>    Time_min Concentration_uM
#> 1         0       15.2808805
#> 2         0       16.0186086
#> 3         0       14.4496542
#> 4        30        8.2960095
#> 5        30        8.7398488
#> 6        30        8.6532021
#> 7        60        5.0325264
#> 8        60        5.1560878
#> 9        60        5.8271061
#> 10      120        0.9598237
#> 11      120        1.3735538
#> 12      120        1.3277293

fit_clearance_from_curve(expData)
#> 85.92071    (5.31e+00): par = (0.01)
#> 6.355533    (1.16e+00): par = (0.01642129)
#> 2.767429    (5.79e-02): par = (0.01858765)
#> 2.758336    (8.26e-04): par = (0.01871273)
#> 2.758334    (1.50e-05): par = (0.01871094)
#> 2.758334    (2.71e-07): par = (0.01871097)

#> Waiting for profiling to be done...
#> Mean_kcat_min-1    2.5%_CI_kcat     95%_CI_kcat 
#>      0.01871097      0.01734030      0.02020731 

```
