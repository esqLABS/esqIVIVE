# IVIVE for clearance

IVIVE for clearance based on different type of values

## Usage

``` r
IVIVE_clearance(
  typeValue,
  units,
  expData,
  typeSystem,
  fu_invitro = 1,
  empirical_scalar = "No",
  tissue = "liver",
  species = "human",
  volMedium_mL = 1,
  REF = 1,
  cProtein_mgml = NULL,
  cCells_Mml = NULL,
  verbose = FALSE
)
```

## Arguments

- typeValue:

  what type of value it is: kcat (directly from in vitro), halfLife or
  invitro_clearance_parameter)

- units:

  this are the units of the value. For kct

  |                         |                       |               |
  |-------------------------|-----------------------|---------------|
  | Cells                   | Subcellular           | Generic       |
  | mL/minutes/millioncells | mL/minutes/mg protein | /minutes      |
  | uL/minutes/millioncells | uL/minutes/mg protein | /hours        |
  | L/minutes/millioncells  | L/minutes/mg protein  | /seconds      |
  | mL/hours/millioncells   | mL/hours/mg protein   | mL/minutes    |
  | uL/hours/millioncells   | uL/hours/mg protein   | uL/minutes    |
  | L/hours/millioncells    | L/hours/mg protein    | mL/seconds    |
  | mL/seconds/millioncells | mL/seconds/mg protein | uL/seconds    |
  | uL/seconds/millioncells | mL/seconds/mg protein | mL/hours      |
  | L/seconds/millioncells  | uL/seconds/mg protein | uL/hours      |
  | mL/minutes/cell         | L/seconds/mg protein  | mL/minutes/kg |
  | uL/minutes/cell         |                       | uL/minutes/kg |
  | mL/hours/cell           |                       | mL/hours/kg   |
  | uL/hours/cell           |                       | uL/hours/kg   |
  | mL/seconds/cell         |                       |               |
  | uL/seconds/cell         |                       |               |

- expData:

  experimental clearance value

- typeSystem:

  "cells", "microsomes" or "cytosolF" (cytosolic fraction). Cells are
  scaled with the cells per gram organ (CellsGO), microsomes with the
  microsomal protein per gram organ (MicProtGO) and the cytosolic
  fraction with the cytosolic protein per gram organ (CytosProtGO)

- fu_invitro:

  Fraction unbound in the in vitro hepatic system. if not known code
  will calculate it

- empirical_scalar:

  this is an option to include an extra empirical correction factor.
  Currently we are considering the scale factor of Wood 2017

- tissue:

  tissue of interest, since scaling factors are dependent on the tissue,
  will default to liver

- species:

  values can human and rat for now, defaulting to human

- volMedium_mL:

  volume of medium in the well (in mL)

- REF:

  relative expression or activity factor

- cProtein_mgml:

  concentration of subcellular (microsomal or cytosolic) protein (in
  mg/mL)

- cCells_Mml:

  concentration of cells used (in million cells/mL)

- verbose:

  if TRUE, print the inputs and the resulting clearance

## Value

Specific clearance parameter (/min) to plug in PK-Sim

## Examples

``` r
# example cells (e.g. hepatocytes)
IVIVE_clearance(typeValue="invitro_clearance_parameter",typeSystem="cells",species="human",
units="mL/minutes/millioncells",expData=18.27,fu_invitro=0.5,cCells_Mml=0.5,empirical_scalar="No")
#> ClspePermin 
#>     6489.94 

# if you dont specify some of the parameters they will be the default (example fu_in vitro=1)
IVIVE_clearance(typeValue="invitro_clearance_parameter",typeSystem="cells",
units="mL/minutes/millioncells",expData=18.27,cCells_Mml=0.5,verbose=TRUE)
#> --- IVIVE_clearance ---
#> Inputs:
#>   typeValue = invitro_clearance_parameter
#>   units = mL/minutes/millioncells
#>   expData = 18.27
#>   typeSystem = cells
#>   fu_invitro = 1
#>   empirical_scalar = No
#>   tissue = liver
#>   species = human
#> Result:
#> Error in .print_ivive_result("IVIVE_clearance", inputs = list(typeValue = typeValue,     units = units, expData = expData, typeSystem = typeSystem,     fu_invitro = fu_invitro, empirical_scalar = empirical_scalar,     tissue = tissue, species = species)): argument "result" is missing, with no default


# example microsomes
IVIVE_clearance(typeValue="invitro_clearance_parameter",typeSystem="microsomes",
                units="L/minutes/mg protein",expData=18.27,fu_invitro=0.5,cProtein_mgml=0.5,
                volMedium_mL=0.5,empirical_scalar="No")
#> ClspePermin 
#>     1963343 

# example cytosolic fraction
IVIVE_clearance(typeValue="invitro_clearance_parameter",typeSystem="cytosolF",
                units="mL/minutes/mg protein",expData=0.05,fu_invitro=0.8,cProtein_mgml=1)
#> ClspePermin 
#>    4.664179 
```
