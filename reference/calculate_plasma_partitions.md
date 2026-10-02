# Calculate the partition coefficients to plasma components

Predicts the partition coefficients of a compound to the albumin, the
globulins and the membrane lipids (such as those of lipoproteins) of
plasma, from its lipophilicity or from its PP-LFER descriptors. Pass the
results to
[`calculate_fu_plasma()`](https://esqlabs.github.io/ESQivive/reference/calculate_fu_plasma.md)
to predict the fraction unbound in plasma.

## Usage

``` r
calculate_plasma_partitions(
  method,
  lipophilicity,
  ionization,
  pka,
  lfer_e = NULL,
  lfer_b = NULL,
  lfer_a = NULL,
  lfer_s = NULL,
  lfer_v = NULL,
  verbose = FALSE
)
```

## Arguments

- method:

  Prediction method: `"logp"` (regressions on lipophilicity) or
  `"pplfer"` (poly-parameter linear free energy relationships).

- lipophilicity:

  Lipophilicity of the compound as logP (log units).

- ionization:

  Ionization class of up to two ionizable groups, as a vector of length
  2 with `"acid"`, `"base"` or `"neutral"`, for example
  `c("acid", "neutral")`.

- pka:

  pKa values of the two ionizable groups, a vector of length 2.

- lfer_e, lfer_b, lfer_a, lfer_s, lfer_v:

  Abraham solute descriptors E, B, A, S and V (the McGowan volume).
  Needed when `method = "pplfer"`.

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

A named list: `partition_albumin` (L/kg), `partition_globulin` (L/kg)
and `partition_membrane_lipids` (L/L).

## Details

For neutral compounds with logP above 4, use `method = "logp"`. For
acidic phenols, carboxylic acids, pyridines and amines you can use
`method = "pplfer"`. Descriptors can be obtained, for example, from the
fup calculator (<https://drumap.nibiohn.go.jp/fup/>).

How ionization is taken into account by the PP-LFER method is still
under review.

## Examples

``` r
calculate_plasma_partitions(
  method = "logp",
  lipophilicity = 2,
  ionization = c("acid", "neutral"),
  pka = c(3, 0)
)
#> $partition_albumin
#> [1] 0.1851041
#> 
#> $partition_globulin
#> [1] 0.3490001
#> 
#> $partition_membrane_lipids
#> [1] 1.000183
#> 

calculate_plasma_partitions(
  method = "pplfer",
  lipophilicity = 2,
  ionization = c("acid", "neutral"),
  pka = c(3, 0),
  lfer_e = 1,
  lfer_b = 0,
  lfer_a = 1.5,
  lfer_s = 0.8,
  lfer_v = 2
)
#> $partition_albumin
#> [1] 16.91123
#> 
#> $partition_globulin
#> [1] 4.852
#> 
#> $partition_membrane_lipids
#> [1] 11324004
#> 
```
