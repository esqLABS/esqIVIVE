# Get measured fraction unbound from the Krumpholz et al. database

Look up measured fraction unbound values (in microsomes, hepatocytes,
plasma or recombinant CYP systems) for one or more compounds in the
Krumpholz et al. dataset shipped in
`inst/extdata/Krumpholz_et_al_fu_dataset.xlsx`. By default all values
measured in the same condition (compound, species and protein or cell
concentration) are averaged. Matching of compound names ignores case and
leading/trailing spaces. If a compound is not in the database a message
is given, with similar compound names when there are any.

## Usage

``` r
get_fu_krumpholz(
  compound,
  system = "microsomes",
  species = NULL,
  average = TRUE,
  verbose = FALSE
)
```

## Arguments

- compound:

  character vector with the compound name(s)

- system:

  in vitro system: "microsomes", "hepatocytes", "plasma" or "recombinant
  CYPs"

- species:

  optional character vector to filter on species (e.g. "human", "rat").
  If NULL all species are returned

- average:

  if TRUE (default), average fu over all records of the same compound,
  species and concentration. If FALSE, return every record with method,
  comments and reference

- verbose:

  if TRUE, print the number of records found per compound

## Value

If `average = TRUE`, a data.frame with one row per condition and columns
compound, species, the concentration of the incubation
(`concentration_mgml` in mg protein/mL for microsomes and recombinant
CYPs, `concentration_Mcellsml` in million cells/mL for hepatocytes, none
for plasma), `cyp` (recombinant CYPs only), `fu` (mean) and `n` (number
of values averaged). Records whose concentration was reported as a range
or not reported are averaged together under an NA concentration.
Duplicate records (same compound, species, concentration, fu and
reference) are counted once.

If `average = FALSE`, a data.frame with one row per record and columns
compound, system, species, method, concentration, concentration_unit,
concentration_reported, fu, fu_sd, fu_range, comments, doi and
reference. For "recombinant CYPs" the columns cyp,
concentration_pmol_cyp_ml and compound_concentration_uM are added.

Compounds that are not found return no rows.

## Examples

``` r
# all verapamil fu_mic in human, averaged per microsomal concentration
get_fu_krumpholz("Verapamil", system = "microsomes", species = "human")
#>    compound species concentration_mgml        fu  n
#> 1 Verapamil   human              0.025 1.0000000  1
#> 2 Verapamil   human              0.250 0.8500000  1
#> 3 Verapamil   human              0.500 0.5808333  6
#> 4 Verapamil   human              0.710 0.5600000  1
#> 5 Verapamil   human              0.760 0.5900000  2
#> 6 Verapamil   human              0.810 0.5250000  1
#> 7 Verapamil   human              1.000 0.5368000 10
#> 8 Verapamil   human              1.490 0.4070000  1
#> 9 Verapamil   human                 NA 0.4300000  1

get_fu_krumpholz(c("midazolam", "Diazepam"), system = "hepatocytes")
#>    compound species concentration_Mcellsml    fu n
#> 1  Diazepam   human                    0.5 0.870 1
#> 2  Diazepam   human                     NA 0.700 2
#> 3  Diazepam     rat                    0.5 0.705 1
#> 4  Diazepam     rat                    1.0 0.840 2
#> 5 Midazolam   human                    0.5 0.710 1
#> 6 Midazolam   human                    1.0 0.360 1
#> 7 Midazolam   human                     NA 0.540 1

# individual records with references
get_fu_krumpholz("Verapamil", system = "microsomes", species = "human", average = FALSE)
#>     compound     system species               method concentration
#> 1  Verapamil microsomes   human equilibrium dialysis         0.760
#> 2  Verapamil microsomes   human equilibrium dialysis            NA
#> 3  Verapamil microsomes   human equilibrium dialysis         1.000
#> 4  Verapamil microsomes   human equilibrium dialysis         0.810
#> 5  Verapamil microsomes   human equilibrium dialysis         1.490
#> 6  Verapamil microsomes   human equilibrium dialysis         0.250
#> 7  Verapamil microsomes   human equilibrium dialysis         0.500
#> 8  Verapamil microsomes   human equilibrium dialysis         0.500
#> 9  Verapamil microsomes   human equilibrium dialysis         1.000
#> 10 Verapamil microsomes   human equilibrium dialysis         1.000
#> 11 Verapamil microsomes   human equilibrium dialysis         0.760
#> 12 Verapamil microsomes   human equilibrium dialysis         1.000
#> 13 Verapamil microsomes   human equilibrium dialysis         1.000
#> 14 Verapamil microsomes   human equilibrium dialysis         1.000
#> 15 Verapamil microsomes   human equilibrium dialysis         1.000
#> 16 Verapamil microsomes   human equilibrium dialysis         1.000
#> 17 Verapamil microsomes   human equilibrium dialysis         0.500
#> 18 Verapamil microsomes   human equilibrium dialysis         0.025
#> 19 Verapamil microsomes   human            HLM-beads         1.000
#> 20 Verapamil microsomes   human            HLM-beads         0.500
#> 21 Verapamil microsomes   human            HLM-beads         0.025
#> 22 Verapamil microsomes   human  ultracentrifugation         0.250
#> 23 Verapamil microsomes   human  ultracentrifugation         1.000
#> 24 Verapamil microsomes   human      ultrafiltration         0.710
#> 25 Verapamil microsomes   human      ultrafiltration         0.500
#> 26 Verapamil microsomes   human      ultrafiltration         1.000
#> 27 Verapamil microsomes   human                 <NA>         0.500
#>    concentration_unit concentration_reported    fu fu_sd fu_range
#> 1       mg protein/mL                   0.76 0.590 0.070     <NA>
#> 2       mg protein/mL                  0.5-1 0.430    NA     <NA>
#> 3       mg protein/mL                      1 0.370    NA     <NA>
#> 4       mg protein/mL                   0.81 0.525    NA     <NA>
#> 5       mg protein/mL                   1.49 0.407    NA     <NA>
#> 6       mg protein/mL                   0.25 0.850    NA     <NA>
#> 7       mg protein/mL                    0.5 0.625 0.076     <NA>
#> 8       mg protein/mL                    0.5 0.630    NA     <NA>
#> 9       mg protein/mL                      1 0.530 0.060     <NA>
#> 10      mg protein/mL                      1 0.470    NA     <NA>
#> 11      mg protein/mL                   0.76 0.590 0.030     <NA>
#> 12      mg protein/mL                      1 0.700 0.040     <NA>
#> 13      mg protein/mL                      1 0.700 0.110     <NA>
#> 14      mg protein/mL                      1 0.730 0.050     <NA>
#> 15      mg protein/mL                      1 0.548    NA     <NA>
#> 16      mg protein/mL                      1 0.430 0.060     <NA>
#> 17      mg protein/mL                    0.5 0.560 0.050     <NA>
#> 18      mg protein/mL  2.5000000000000001E-2 1.000 0.100     <NA>
#> 19      mg protein/mL                      1 0.320 0.320     <NA>
#> 20      mg protein/mL                    0.5 0.590 0.040     <NA>
#> 21      mg protein/mL  2.5000000000000001E-2 1.000 0.100     <NA>
#> 22      mg protein/mL                   0.25 0.850    NA     <NA>
#> 23      mg protein/mL                      1 0.670 0.100     <NA>
#> 24      mg protein/mL                   0.71 0.560    NA     <NA>
#> 25      mg protein/mL                    0.5 0.650    NA     <NA>
#> 26      mg protein/mL                      1 0.600 0.020     <NA>
#> 27      mg protein/mL                    0.5 0.430 0.100     <NA>
#>                                                                             comments
#> 1  final protein concentration of 0.76 mg/mL. final incubation concentration of 1 μM
#> 2                                                                        0.5-1 mg/ml
#> 3                                                                               <NA>
#> 4                                                                               <NA>
#> 5                                                                               <NA>
#> 6                                                                               <NA>
#> 7                                                                               <NA>
#> 8                                                                            Tabl. 2
#> 9                                                                            1 mg/ml
#> 10                                                                         1.0 mg/ml
#> 11                                      1 μM; final protein concentration 0.76 mg/ml
#> 12                                                                   1 mg/ml; 100 μM
#> 13                                                                   1 mg/ml; 500 μM
#> 14                                                                   1 mg/ml; 200 μM
#> 15                                                                              <NA>
#> 16                                                                              <NA>
#> 17                                                                              <NA>
#> 18                                                                              <NA>
#> 19                                                                              <NA>
#> 20                                                                              <NA>
#> 21                                                                              <NA>
#> 22                                                                              <NA>
#> 23                                                                           1 mg/ml
#> 24                                                            HLM protein 0.71 mg/ml
#> 25                                                             HLM protein 0.5 mg/ml
#> 26                                                                           1 mg/ml
#> 27                                                microsomal concentration 0.5 mg/ml
#>                                                                                                                                        doi
#> 1                                                                                                                        10.1002/jps.22124
#> 2                                                                                                                10.1007/s11095-006-9663-4
#> 3                                                                                                               10.1007/s11095-022-03205-1
#> 4                                                                                                               10.1007/s11095-022-03205-1
#> 5                                                                                                              10.1016/j.vascn.2005.10.002
#> 6                                                                                                              10.1016/j.vascn.2010.04.003
#> 7                                                                                                               10.1016/j.xphs.2020.09.012
#> 8                                                                                                            10.1080/00498254.2022.2132426
#> 9                                                                                                                   10.1124/dmd.105.005033
#> 10                                                                                                                  10.1124/dmd.107.018713
#> 11                                                                                                                  10.1124/dmd.107.020131
#> 12                                                                                                                  10.1124/dmd.111.039354
#> 13                                                                                                                  10.1124/dmd.111.039354
#> 14                                                                                                                  10.1124/dmd.111.039354
#> 15                                                                                                                  10.1124/dmd.120.000131
#> 16                                                                                                                  10.1124/dmd.121.000575
#> 17                                                                                                                  10.1124/dmd.121.000575
#> 18                                                                                                                  10.1124/dmd.121.000575
#> 19                                                                                                                  10.1124/dmd.121.000575
#> 20                                                                                                                  10.1124/dmd.121.000575
#> 21                                                                                                                  10.1124/dmd.121.000575
#> 22                                                                                                             10.1016/j.vascn.2010.04.003
#> 23                                                                                                                  10.1124/dmd.105.005033
#> 24                                                                                                                       10.1002/jps.22635
#> 25                                                                                                                       10.1002/jps.22635
#> 26                                                                                                                  10.1124/dmd.105.005033
#> 27 https://www.semanticscholar.org/paper/Prediction-of-human-clearance-of-twenty-nine-drugs-Obach/c1bded353d7f63823e4a4fa21a872f77f34d55d3
#>          reference
#> 1       Zhang 2010
#> 2    Mohutsky 2006
#> 3        Tess 2022
#> 4        Tess 2022
#> 5      Skaggs 2006
#> 6    Deshmukh 2011
#> 7       Jones 2021
#> 8     Gardner 2022
#> 9    Giuliano 2005
#> 10      Gertz 2008
#> 11        Gao 2008
#> 12     McLure 2011
#> 13     McLure 2011
#> 14     McLure 2011
#> 15 Williamson 2020
#> 16       Wang 2021
#> 17       Wang 2021
#> 18       Wang 2021
#> 19       Wang 2021
#> 20       Wang 2021
#> 21       Wang 2021
#> 22   Deshmukh 2011
#> 23   Giuliano 2005
#> 24   Beaumont 2011
#> 25   Beaumont 2011
#> 26   Giuliano 2005
#> 27      Obach 1999

# check if a compound is in the database
nrow(get_fu_krumpholz("Verapamil", system = "plasma")) > 0
#> [1] TRUE
```
