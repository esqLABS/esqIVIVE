# Calculate Fu in plasma based on affinity to different components

Calculate the Fu in plasma based on the affinity to the different
components: albumin, globulin and membrane lipids such as the ones in
lipoproteins

## Usage

``` r
calculate_fu_pls_from_Ks(
  partition_albumin,
  partition_globulin,
  partition_membrane_lipids,
  partition_lipids,
  species,
  verbose = FALSE
)
```

## Arguments

- partition_albumin:

  partition to albumin (in L/kg)

- partition_globulin:

  partition to globulin (in L/kg)

- partition_membrane_lipids:

  partition to membrane lipids (in L/L)

- partition_lipids:

  LogP or LogMA of the compound

- species:

  species to be considered, now there is data for human, rat, dog,
  monkey, rabbit and mouse. Mind that this changes the composition in
  serum but the specific affinities to albumin still need to be used

- verbose:

  if TRUE, print the inputs and the resulting Fu_plasma

## Value

Fu_plasma value

## Examples

``` r
calculate_fu_pls_from_Ks("partition_albumin"=10^4.48,
               "partition_globulin"=10^2.16,
               "partition_membrane_lipids"=10^3.51,
               "partition_lipids"=100,
               "species"="human")
#>   Fu_plasma 
#> 0.000798041 
```
