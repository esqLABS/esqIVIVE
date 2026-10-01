# Poulin algorithm for Fu calculation

The Poulin algorithm for calculating Fu in vitro

## Usage

``` r
calculate_fu_hep_poulin(
  ionization,
  pKa,
  blood_plasma = NULL,
  fraction_unbound = NULL,
  concentration_cell_neutral_lipids,
  cCellAPL,
  log_lipophilicity,
  verbose = FALSE
)
```

## Arguments

- ionization:

  Vector of length 2 with ionization class, acid, neutral and base, if
  not input then it is c(0,0)

- pKa:

  vector of length of 2 with pKa of the compound

- blood_plasma:

  blood plasma ratio

- fraction_unbound:

  fraction unbound in plasma

- concentration_cell_neutral_lipids:

  neutral lipid concentration

- cCellAPL:

  acidic phospholipid concentration in the in vitro system (as fraction
  of medium volume); only used for strong bases (base ionization class
  with a pKa above 7)

- log_lipophilicity:

  LogP or LogMA of the compound

- verbose:

  if TRUE, print the inputs and the resulting fu_invitro

## Value

fuInvitro

## Examples

``` r
calculate_fu_hep_poulin(ionization=c("base",0),pKa=c(8,0),blood_plasma=1,fraction_unbound=0.2,
                        concentration_cell_neutral_lipids=0.03,cCellAPL=0.01,log_lipophilicity=3)
#> [1] 0.02036919
calculate_fu_hep_poulin(ionization=c("neutral",0),pKa=c(0,0),
                        concentration_cell_neutral_lipids=0.03,log_lipophilicity=3)
#> [1] 0.03225806
```
