
<!-- README.md is generated from README.Rmd. Please edit that file -->

# esqIVIVE

<!-- badges: start -->
<!-- badges: end -->

The goal of esqIVIVE is to perform extrapolation of in vitro ADME parameters 
and in vitro effect concentrations to in vivo. 
The herein extrapolations have been developed focusing on the integration with OSP tools.

In the near future the codes will be transformed in a package . 
As of now the public can use the code as is, examples are provided.
Currently there are available codes to calculate:
- the fraction unbound in in vitro hepatic models (microsomes and hepatocytes)
- the calculation or extrapolation of in vitro hepatic clearance : https://github.com/esqLABS/esqIVIVE/blob/main/R/clearance_IVIVE.R
- the calculation or extrapolation of in vitro hepatic metabolic Vmax and Km
- the calculation of fraction unbound in plasma using protein-specific affinity constant


## Installation

You can install the development version of esqIVIVE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("esqLABS/esqIVIVE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(esqIVIVE)
## basic example code
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

To run packages tests, execute

``` r
devtools::test()
```

### Documentation Generation

``` r
pkgdown:::build_site_external()
```
