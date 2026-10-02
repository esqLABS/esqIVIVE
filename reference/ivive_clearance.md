# Scale an in vitro clearance to in vivo

Converts an in vitro half-life, depletion rate constant or intrinsic
clearance into the specific clearance (1/min) that PK-Sim expects. The
value is corrected for binding in the incubation and scaled to the whole
tissue with the physiological scaling factors of the species.

## Usage

``` r
ivive_clearance(
  value_type,
  value,
  unit,
  system,
  fu_in_vitro = 1,
  concentration_microsomes = NULL,
  concentration_cells = NULL,
  concentration_cytosol = NULL,
  volume_medium = 1,
  empirical_correction = FALSE,
  tissue = "liver",
  species = "human",
  relative_expression_factor = 1,
  verbose = FALSE
)
```

## Arguments

- value_type:

  Type of `value`:

  - `"half_life"`: in vitro half-life of substrate depletion;

  - `"rate_constant"`: depletion rate constant, for example from
    [`fit_depletion_curve()`](https://esqlabs.github.io/ESQivive/reference/fit_depletion_curve.md);

  - `"intrinsic_clearance"`: in vitro intrinsic clearance.

- value:

  Measured value, in `unit`.

- unit:

  Unit of `value`. The allowed units depend on `value_type`:

  - `"half_life"`: `"minutes"`, `"hours"`, `"seconds"`.

  - `"rate_constant"`: `"/minutes"`, `"/hours"`, `"/seconds"`.

  - `"intrinsic_clearance"` per million cells, per cell or per mg
    protein: a volume (`mL`, `uL`, `L`), a time (`minutes`, `hours`,
    `seconds`) and the amount, for example `"mL/minutes/millioncells"`,
    `"uL/hours/cell"` or `"L/minutes/mg protein"`. `L` is not available
    per cell.

  - `"intrinsic_clearance"` per incubation: `"mL/minutes"`,
    `"uL/minutes"`, `"mL/seconds"`, `"uL/seconds"`, `"mL/hours"`,
    `"uL/hours"`.

  - `"intrinsic_clearance"` per kg body weight: `"mL/minutes/kg"`,
    `"uL/minutes/kg"`, `"mL/hours/kg"`, `"uL/hours/kg"`.

- system:

  Incubation system: `"microsomes"`, `"cells"` (for example hepatocytes)
  or `"cytosol"` (cytosolic fraction). Microsomes are scaled with the
  microsomal protein per gram tissue, cells with the cells per gram
  tissue and the cytosol with the cytosolic protein per gram tissue.

- fu_in_vitro:

  Fraction unbound in the incubation, for example from
  [`calculate_fu_in_vitro()`](https://esqlabs.github.io/ESQivive/reference/calculate_fu_in_vitro.md).
  Defaults to 1 (no binding).

- concentration_microsomes:

  Microsomal protein concentration (mg/mL). Needed for microsomes when
  `value_type` is `"half_life"` or `"rate_constant"`, or when `unit` is
  per incubation.

- concentration_cells:

  Cell concentration (million cells/mL). Needed for cells in the same
  cases.

- concentration_cytosol:

  Cytosolic protein concentration (mg/mL). Needed for the cytosol in the
  same cases.

- volume_medium:

  Volume of medium in the incubation (mL), used for the per incubation
  units. Defaults to 1.

- empirical_correction:

  If `TRUE`, apply the empirical correction factors of Wood et al.
  (2017), which correct the tendency of in vitro data to overpredict
  slow and underpredict fast clearances. Available for human and rat,
  with microsomes or cells.

- tissue:

  Tissue whose scaling factors are used. Defaults to `"liver"`.

- species:

  Species whose scaling factors are used: `"human"`, `"rat"` or `"dog"`.
  Defaults to `"human"`.

- relative_expression_factor:

  Relative expression or activity factor of the enzyme in vivo compared
  with the incubation. Defaults to 1.

- verbose:

  If `TRUE`, print the inputs.

## Value

The specific clearance (1/min), a single number.

## Examples

``` r
# cells, for example hepatocytes
ivive_clearance(
  value_type = "intrinsic_clearance",
  value = 18.27,
  unit = "mL/minutes/millioncells",
  system = "cells",
  fu_in_vitro = 0.5,
  concentration_cells = 0.5
)
#> [1] 6489.94

# microsomes, from the in vitro half-life
ivive_clearance(
  value_type = "half_life",
  value = 3.9,
  unit = "minutes",
  system = "microsomes",
  fu_in_vitro = 0.4,
  concentration_microsomes = 1
)
#> [1] 23.86912

# cytosolic fraction
ivive_clearance(
  value_type = "intrinsic_clearance",
  value = 0.05,
  unit = "mL/minutes/mg protein",
  system = "cytosol",
  fu_in_vitro = 0.8,
  concentration_cytosol = 1
)
#> [1] 4.664179
```
