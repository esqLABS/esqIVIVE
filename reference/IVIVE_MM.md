# Title

function that scales Vmax and correct km for fraction unbound

## Usage

``` r
IVIVE_MM(
  typeSystem,
  fu_invitro = 1,
  vmax,
  km_micromolar,
  tissue = "liver",
  species = "human",
  REF = 1,
  verbose = FALSE
)
```

## Arguments

- typeSystem:

  if hepatocytes or microsomes

- fu_invitro:

  value of fractionunbound in vitro, the default is 1

- vmax:

  as umol/min/million hepatocytes or umol/min/mg microsomal protein

- km_micromolar:

  Km of the enzyme reaction, in uM

- tissue:

  liver, brain, lung, kidney, gonads and gut, default is liver

- species:

  human, rat or dog, default human

- REF:

  relative expression or activity factor, default is 1. To use this
  option the reference concentration of of the enzyme of interest in
  pksim needs to be 1 uM

- verbose:

  if TRUE, print the inputs and the resulting Vmax/Km

## Value

Vmax in umol/min/L and Km_unb in uM

## Examples

``` r
IVIVE_MM (typeSystem="hepatocytes",vmax=2,km_micromolar=1,tissue="liver",species="human",REF=1)
#> $vmax_umol_minL
#> [1] 355223.9
#> 
#> $Km_unb_uM
#> [1] 1
#> 
IVIVE_MM (typeSystem="microsomes",fu_invitro=0.2,vmax=2,km_micromolar=1)
#> $vmax_umol_minL
#> [1] 107462.7
#> 
#> $Km_unb_uM
#> [1] 0.2
#> 
```
