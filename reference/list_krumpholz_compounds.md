# List compounds in the Krumpholz et al. database

List compounds in the Krumpholz et al. database

## Usage

``` r
list_krumpholz_compounds(system = "microsomes")
```

## Arguments

- system:

  in vitro system: "microsomes", "hepatocytes", "plasma" or "recombinant
  CYPs"

## Value

sorted character vector with the compound names available for that
system

## Examples

``` r
head(list_krumpholz_compounds("hepatocytes"))
#> [1] "Acetaminophen" "Albendazole"   "Aldosterone"   "Alprazolam"   
#> [5] "Amiodarone"    "Amitriptyline"
```
