# Austin algorithm for hepatocytes Fu calculation

The Austin algorithm for calculating Fu in vitro for hepatocytes

## Usage

``` r
calculate_fu_hep_austin(
  ionization,
  pKa,
  log_lipophilicity,
  conc_cell_millionml,
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

- conc_cell_millionml:

  concentration of hepatocytes (in million cells/mL)

- verbose:

  if TRUE, print the inputs and the resulting fu_invitro

## Value

fuInvitro

## Examples

``` r
calculate_fu_hep_austin(ionization=c("base",0),pKa=c(3,0),log_lipophilicity=3,
                        conc_cell_millionml=0.5)
#> [1] 0.7516837
```
