# Halifax algorithm for Fu calculation

The Halifax algorithm for calculating Fu in vitro

## Usage

``` r
calculate_fu_mic_halifax(
  ionization,
  pKa,
  log_lipophilicity,
  conc_mic_mgml,
  verbose = FALSE
)
```

## Arguments

- ionization:

  Vector of length 2 with ionization class, acid, neutral and base, if
  not input then it is c(0,0)

- pKa:

  vector of length of 2 with pKa of the compound

- log_lipophilicity:

  LogP or LogMA of the compound

- conc_mic_mgml:

  concentration of microsomes mg/mL

- verbose:

  if TRUE, print the inputs and the resulting fu_invitro

## Value

fuInvitro

## Examples

``` r
calculate_fu_mic_halifax(c("base",0),c(3,0),3,1)
#> [1] 0.6542596
```
