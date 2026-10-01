# Predict affinity constant to plasma components based on QSARs

Collection of QSAR to obtain affinity to the different components in
serum, membrane lipids (memlip) albumin (alb) and globulin (glob)

## Usage

``` r
predict_plasma_affinities(
  QSAR,
  logP,
  pKa,
  ionization,
  LFER_E = NULL,
  LFER_B = NULL,
  LFER_A = NULL,
  LFER_S = NULL,
  LFER_V = NULL,
  verbose = FALSE
)
```

## Arguments

- QSAR:

  type of QSAR, it can be logP based or PPLFER based (still sorting the
  ionization)

- logP:

  is the lipophilicity as given by logKow

- pKa:

  is a vector of length of 2

- ionization:

  is a vector of length of 2 which should indicate if chemical is
  neutral, acid basic. There are spots for ionization in case chemical
  is zwitterion

- LFER_E:

  LFER E parameter

- LFER_B:

  LFER B parameter

- LFER_A:

  LFER A parameter

- LFER_S:

  LFER S parameter

- LFER_V:

  abraham volume

- verbose:

  if TRUE, print the inputs and the resulting partition coefficients

## Value

partition_membrane_lipids (in L/L), partition_albumin (in L/kg) and
partition_globulin (in L/kg)

## Details

For neutral chemicals with logP \>4 use the logP QSAR For acidic
phenols, carboxylic acids, pyridine and amines you can use the PPLFER.
fup calculator (https://drumap.nibiohn.go.jp/fup/). To Do: PP-LFER
QSARs.-need to check how ionization is considered make documentation

## Examples

``` r

predict_plasma_affinities(QSAR="logP", logP=2, pKa=c(3,0), ionization=c("acid",0))
#> partition_membrane_lipids.ion_factor_plasma 
#>                                   1.0001833 
#>         partition_albumin.ion_factor_plasma 
#>                                   0.1851041 
#>                          partition_globulin 
#>                                   0.3490001 

predict_plasma_affinities(QSAR="PPLFER", logP=2, pKa=c(3,0), ionization=c("acid",0),
                          LFER_E=1, LFER_B=0, LFER_A=1.5, LFER_S=0.8, LFER_V=2)
#> partition_membrane_lipids         partition_albumin        partition_globulin 
#>              1.132400e+07              1.691123e+01              4.852000e+00 
```
