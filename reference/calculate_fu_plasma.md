# Calculate the fraction unbound in plasma from partition coefficients

Calculates the fraction unbound in plasma from the partition
coefficients of the compound to albumin, globulins, membrane lipids
(such as those of lipoproteins) and neutral lipids, and the plasma
composition of the species. The first three partition coefficients can
be predicted with
[`calculate_plasma_partitions()`](https://esqlabs.github.io/ESQivive/reference/calculate_plasma_partitions.md).

## Usage

``` r
calculate_fu_plasma(
  partition_albumin,
  partition_globulin,
  partition_membrane_lipids,
  partition_neutral_lipids,
  species,
  verbose = FALSE
)
```

## Arguments

- partition_albumin:

  Partition coefficient to albumin (L/kg).

- partition_globulin:

  Partition coefficient to globulins (L/kg).

- partition_membrane_lipids:

  Partition coefficient to membrane lipids (L/L).

- partition_neutral_lipids:

  Partition coefficient to neutral lipids (L/L). The value is used as
  given, not as a log value.

- species:

  Species whose plasma composition is used: `"human"`, `"rat"`, `"dog"`,
  `"monkey"`, `"rabbit"` or `"mouse"`. Only the plasma composition
  changes with the species: the partition coefficients must be those of
  that species.

- verbose:

  If `TRUE`, print the inputs and the result, with the fractions bound
  to albumin, globulins and lipids.

## Value

The fraction unbound in plasma, a single number.

## Examples

``` r
calculate_fu_plasma(
  partition_albumin = 10^4.48,
  partition_globulin = 10^2.16,
  partition_membrane_lipids = 10^3.51,
  partition_neutral_lipids = 100,
  species = "human"
)
#> [1] 0.000798041
```
