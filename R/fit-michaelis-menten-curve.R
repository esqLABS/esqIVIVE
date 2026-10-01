#' Fit a Michaelis-Menten curve
#'
#' @description
#' Fits the Michaelis-Menten equation to reaction velocities measured at
#' several substrate concentrations and returns Km and Vmax. A plot of the
#' data and the fitted curve is drawn so you can judge the fit.
#'
#' To scale the results to in vivo, pass them to [ivive_michaelis_menten()].
#'
#' @param data A data frame with the substrate concentration (uM) in the first
#'   column and the velocity in the second column.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @returns A data frame with two rows, `parameter = "km"` (uM) and
#'   `parameter = "vmax"` (in the velocity unit of `data`), and the columns
#'   `estimate`, `lower` and `upper`, the bounds of the 95% confidence
#'   interval. A warning is given when the fit is poor (R-squared below 0.8).
#' @export
#'
#' @examples
#' mm_curve <- read.csv(
#'   system.file("extdata", "michaelis_menten_curve.csv", package = "ESQivive")
#' )
#' fit_michaelis_menten_curve(mm_curve)
fit_michaelis_menten_curve <- function(data, verbose = FALSE) {
  experimental_conc_velocity <- data[, 1:2]
  colnames(experimental_conc_velocity) <- c("Concentration", "Velocity")

  #fit model
  fitmm <- stats::nls(
    Velocity ~ Vmax * Concentration / (Km + Concentration),
    data = experimental_conc_velocity,
    start = list(
      Vmax = max(experimental_conc_velocity$Velocity),
      Km = mean(experimental_conc_velocity$Concentration)
    ),
    trace = FALSE
  )

  rss <- sum(stats::residuals(fitmm)^2)
  tss <- sum(
    (experimental_conc_velocity$Velocity -
      mean(experimental_conc_velocity$Velocity))^2
  )
  r2 <- round(1 - (rss / tss), digits = 3)

  if (r2 < 0.8) {
    cli::cli_warn(
      "Poor fit of the Michaelis-Menten curve (R-squared {r2}): the data may
       not follow Michaelis-Menten kinetics."
    )
  }

  fit_95conf <- stats::confint(fitmm)

  #check if fitting is good
  mm_fuction <- function(Concentration) {
    stats::coefficients(fitmm)[["Vmax"]] *
      Concentration /
      (stats::coefficients(fitmm)[["Km"]] + Concentration)
  }

  plot_diagnosis <- ggplot2::ggplot(
    experimental_conc_velocity,
    ggplot2::aes(x = .data$Concentration, y = .data$Velocity)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(x = colnames(data)[1], y = colnames(data)[2]) +
    ggplot2::stat_function(fun = mm_fuction, colour = "blue") +
    ggplot2::annotate(
      "text",
      x = max(experimental_conc_velocity$Concentration) * 0.8,
      y = max(experimental_conc_velocity$Velocity) * 0.8,
      label = paste("R\u00b2=", r2),
      size = 5
    )

  #Add fit values in dataframe for calculations
  result <- data.frame(
    parameter = c("km", "vmax"),
    estimate = unname(stats::coefficients(fitmm)[c("Km", "Vmax")]),
    lower = unname(fit_95conf[c("Km", "Vmax"), 1]),
    upper = unname(fit_95conf[c("Km", "Vmax"), 2])
  )

  print(plot_diagnosis)

  if (verbose) {
    .print_ivive_result(
      "fit_michaelis_menten_curve",
      inputs = list(data = data),
      result = result
    )
  }

  result
}
