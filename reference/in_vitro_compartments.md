# in_vitro_compartments

Generates a list of values describing a liver in vitro compartment based
on hepatocytes or microsomes This function is used inside the Fraction
unbound function but can also be used for general virtual hepatocyte
systems

## Usage

``` r
in_vitro_compartments(
  typeSystem,
  FBS_fraction,
  microplateType,
  volMedium_mL,
  cCells_Mml = NULL,
  cMicro_mgml = NULL,
  verbose = FALSE
)
```

## Arguments

- typeSystem:

  if system is hepatocytes or microsomes

- FBS_fraction:

  fraction of serum concentration, values can only go from 0-1

- microplateType:

  number of wells in the microplate

- volMedium_mL:

  volume of medium in the well (in mL)

- cCells_Mml:

  cells concentration (in million cells/mL)

- cMicro_mgml:

  concentration of microsome protein (in mg/mL)

- verbose:

  if TRUE, print the inputs and the resulting compartment list

## Value

a list of values representing the different in vitro compartments,
concentrations are given as fraction of volume

## Examples

``` r
in_vitro_compartments("hepatocytes", FBS_fraction=0.05, microplateType = 96,
                      volMedium_mL = 0.15, cCells_Mml = 0.1)
#> $cCellNL_vvmedium
#> [1] 1.1303e-05
#> 
#> $cCellNPL_vvmedium
#> [1] 8.4074e-06
#> 
#> $cCellAPL_vvmedium
#> [1] 2.2352e-06
#> 
#> $cCellPro_vvmedium
#> [1] 5.08e-05
#> 
#> $cMediumNL_vvmedium
#> [1] 7.85e-05
#> 
#> $cMediumNPL_vvmedium
#> [1] 1.5e-05
#> 
#> $cMediumPro_vvmedium
#> [1] 0.002
#> 
#> $saPlasticVolMedium_m2L
#> [1] 0.6060606
#> 
#> $volAir_L
#> [1] 0.000242
#> 
in_vitro_compartments("microsomes", FBS_fraction=0, microplateType = 24,
                      volMedium_mL = 0.5, cMicro_mgml = 1)
#> $cCellNL_vvmedium
#> [1] 0.0002611111
#> 
#> $cCellNPL_vvmedium
#> [1] 0.0007261556
#> 
#> $cCellAPL_vvmedium
#> [1] 0.0001594
#> 
#> $cCellPro_vvmedium
#> [1] 0.0007407407
#> 
#> $cMediumNL_vvmedium
#> [1] 0
#> 
#> $cMediumNPL_vvmedium
#> [1] 0
#> 
#> $cMediumPro_vvmedium
#> [1] 0
#> 
#> $saPlasticVolMedium_m2L
#> [1] 0
#> 
#> $volAir_L
#> [1] 0.00297
#> 
```
