---
editor_options: 
  markdown: 
    wrap: 72
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# esqIVIVE

<!-- badges: start -->

<!-- badges: end -->

The goal of esqIVIVE is to perform extrapolation of in vitro ADME
parameters and derive ADME parameters to input for PBK models.

The functions in this package have been developed focusing on the
integration with OSP tools.

Currently there are available codes to calculate,

[fraction unbound in microsomes]{.underline} :

-   calculate_fu_mic_austin()

-   calculate_fu_hep_halifax()

-   calculate_fu_mic_turner()

[fraction unbound in hepatocytes]{.underline}:

-   calculate_fu_hep_austin()

-   calculate_fu_hep_kilford()

-   calculate_fu_hep_poulin()

[derive metabolism parameters from experimental curves:]{.underline}

-   fit_clearance_from_curve()

-   get_MM()

[perform scaling for clearance:]{.underline}

-   IVIVE_clearance()
-   IVIVE_MM()

[calculate fu_plasma related parameters:]{.underline}

-   calculate_fu_pls_from_Ks()

-   predict_plasma_affinities()

-   correct_fu_pls_pearce()

[perform IVIVE to derive Pint:]{.underline}

-   pint_caco2_empir()

-   pint_peff_empir()

Examples of how to use the functions are provided for each.

## Installation

You can install the development version of esqIVIVE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("esqLABS/esqIVIVE")
```

## Contribute

### Coding Standards

Contributors should comply with the [Open Systems Pharmacology Coding
Standards for
R](https://github.com/Open-Systems-Pharmacology/developer-docs/blob/main/ospsuite-r-specifics/CODING_STANDARDS_R.md)

### Development Environment

To install all the dependencies required for development, run:

``` r
renv::install()
```

### Testing

To run packages tests, execute ( dont do for now)

``` r
devtools::test()
```

### Documentation Generation

``` r
pkgdown:::build_site_external()
```
