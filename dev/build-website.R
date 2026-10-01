# Run from the package root: Rscript dev/build-website.R
# Quarto and the ospsuite/PK-Sim system prerequisites must already be installed.
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}

pak::pkg_install("pkgdown", ask = FALSE)
pak::local_install(
  dependencies = c("Depends", "Imports", "LinkingTo", "Suggests", "Config/Needs/website"),
  upgrade = FALSE,
  ask = FALSE
)

pkgdown::build_site()
