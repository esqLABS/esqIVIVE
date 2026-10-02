# Calculate the fraction unbound in an in vitro incubation

Predicts the fraction unbound of a compound in an incubation with
microsomes or cells, such as hepatocytes. Two kinds of method are
available:

- literature regressions on lipophilicity (`austin`, `hallifax`,
  `turner`, `kilford`, `poulin`, and their average `all_literature`),
  which need only the compound properties and the microsome or cell
  concentration;

- partition models (`poulin_theil`, `berezhkovskiy`, `pksim_standard`,
  `rodgers_rowland`, `schmitt`), which describe binding to the lipids
  and proteins of the cells or microsomes, to serum in the medium and to
  the plastic of the well. They also need the serum fraction, the
  microplate type and the medium volume.

## Usage

``` r
calculate_fu_in_vitro(
  method,
  system,
  lipophilicity,
  ionization,
  pka = NULL,
  concentration_microsomes = NULL,
  concentration_cells = NULL,
  fbs_fraction = NULL,
  microplate_type = NULL,
  volume_medium = NULL,
  henry_law_constant = NULL,
  fu_plasma = NULL,
  blood_plasma_ratio = NULL,
  verbose = FALSE
)
```

## Arguments

- method:

  Prediction method, one of:

  |  |  |  |
  |----|----|----|
  | `method` | System | Description |
  | `"austin"` | both | Austin et al. (2002) regression |
  | `"hallifax"` | microsomes | Hallifax and Houston (2006) regression |
  | `"turner"` | microsomes | Turner regression, with separate equations for acids, bases and neutral compounds |
  | `"kilford"` | cells | Kilford et al. (2008) regression |
  | `"poulin"` | both | Poulin regression on the neutral lipid content, with acidic phospholipid binding for strong bases |
  | `"all_literature"` | both | Average of the regressions available for the system |
  | `"poulin_theil"`, `"poulin_theil_fu"` | both | Poulin and Theil partition model |
  | `"berezhkovskiy"`, `"berezhkovskiy_fu"` | both | Berezhkovskiy partition model |
  | `"pksim_standard"`, `"pksim_standard_fu"` | both | PK-Sim Standard partition model |
  | `"rodgers_rowland_fu"` | both | Rodgers and Rowland partition model, strong bases only |
  | `"schmitt"`, `"schmitt_fu"` | both | Schmitt partition model |

  The partition models without the `_fu` suffix predict binding to serum
  in the medium from the serum lipid and protein content. The `_fu`
  versions use the measured `fu_plasma` instead.

- system:

  Incubation system, `"microsomes"` or `"cells"` (for example
  hepatocytes).

- lipophilicity:

  Lipophilicity of the compound (log units), as logP or log membrane
  affinity.

- ionization:

  Ionization class of up to two ionizable groups, as a vector of length
  2 with `"acid"`, `"base"` or `"neutral"`, for example
  `c("base", "neutral")`.

- pka:

  pKa values of the two ionizable groups, a vector of length 2. Defaults
  to `c(0, 0)`.

- concentration_microsomes:

  Microsomal protein concentration (mg/mL). Needed when
  `system = "microsomes"`.

- concentration_cells:

  Cell concentration (million cells/mL). Needed when `system = "cells"`.

- fbs_fraction:

  Fraction of fetal bovine serum in the medium (0 to 1). Needed by the
  partition models.

- microplate_type:

  Number of wells of the microplate: 96, 48, 24 or 12. Needed by the
  partition models.

- volume_medium:

  Volume of medium in the well (mL). Needed by the partition models.

- henry_law_constant:

  Henry's law constant (atm m3/mol), used to warn when the compound
  probably evaporates from the well. Defaults to 1e-6.

- fu_plasma:

  Fraction unbound in plasma. Needed by the `_fu` methods, and by
  `"poulin"` and `"all_literature"` for strong bases.

- blood_plasma_ratio:

  Blood to plasma concentration ratio. Needed by `"rodgers_rowland_fu"`,
  and by `"poulin"` and `"all_literature"` for strong bases.

- verbose:

  If `TRUE`, print the inputs and the result.

## Value

The fraction unbound in the incubation, a single number between 0 and 1.
A warning is given when more than 5% of the compound is predicted to be
in the air of the well; this check needs `microplate_type` and
`volume_medium`.

## Details

A strong base is a compound whose first ionization class is `"base"`
with a pKa above 7.

`"berezhkovskiy"` and `"berezhkovskiy_fu"` currently give the same
results as `"poulin_theil"` and `"poulin_theil_fu"`.

## Examples

``` r
calculate_fu_in_vitro(
  method = "austin",
  system = "microsomes",
  lipophilicity = 3,
  ionization = c("base", "neutral"),
  pka = c(8, 0),
  concentration_microsomes = 1
)
#> [1] 0.5689242

calculate_fu_in_vitro(
  method = "pksim_standard",
  system = "cells",
  lipophilicity = 3,
  ionization = c("acid", "neutral"),
  pka = c(6, 0),
  concentration_cells = 2,
  fbs_fraction = 0,
  microplate_type = 96,
  volume_medium = 0.22
)
#> [1] 0.5630088

calculate_fu_in_vitro(
  method = "poulin_theil_fu",
  system = "cells",
  lipophilicity = 3,
  ionization = c("acid", "neutral"),
  pka = c(6, 0),
  concentration_cells = 2,
  fbs_fraction = 0.05,
  microplate_type = 96,
  volume_medium = 0.22,
  fu_plasma = 0.01
)
#> [1] 0.1606004
```
