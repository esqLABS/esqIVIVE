#' Fit a substrate depletion curve
#'
#' @description
#' Fits a one-phase exponential decay to the concentrations measured in a
#' substrate depletion experiment and returns the depletion rate constant.
#' The starting concentration is the mean of the concentrations at time 0.
#' A plot of the data and the fitted curve is drawn so you can judge the fit.
#'
#' The rate constant is not yet a clearance for PK-Sim: pass it to
#' [ivive_clearance()] with `value_type = "rate_constant"`.
#'
#' @param data A data frame with the time (min) in the first column and the
#'   concentration (uM) in the second column. Include the time 0 samples.
#' @param verbose If `TRUE`, print the inputs and the result.
#'
#' @returns A data frame with one row, `parameter = "rate_constant"` (1/min),
#'   and the columns `estimate`, `lower` and `upper`, the bounds of the 95%
#'   confidence interval. A warning is given when the fit is poor (R-squared
#'   below 0.8) or when the concentration falls by less than 80%.
#' @export
#'
#' @examples
#' depletion <- read.csv(system.file("extdata", "clearance.csv", package = "ESQivive"))
#' head(depletion)
#'
#' fit_depletion_curve(depletion)
fit_depletion_curve <- function(data, verbose = FALSE) {
  #Load the depletion curve
  clear_curve_xy <- data[, 1:2]
  colnames(clear_curve_xy) <- c("x", "y")

  #find the starting concentration
  y0 <- mean(clear_curve_xy$y[clear_curve_xy$x == 0])

  #create clearance model
  Kcat_function <- function(x, clearance_rate_constant) {
    y0 * exp(-clearance_rate_constant * x)
  }

  #fit model
  fitKcat <- stats::nls(
    y ~ Kcat_function(x, clearance_rate_constant),
    data = clear_curve_xy,
    start = list(clearance_rate_constant = 0.01),
    trace = TRUE
  )

  rss <- sum(stats::residuals(fitKcat)^2)
  tss <- sum((clear_curve_xy$y - mean(clear_curve_xy$y))^2)
  r2 <- round(1 - (rss / tss), digits = 3)

  #Plot for evaluating if fit is reasonable
  diag_plot <- ggplot2::ggplot(
    clear_curve_xy,
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "fit curve",
      x = colnames(data)[1],
      y = colnames(data)[2]
    ) +
    ggplot2::stat_function(
      fun = function(x) {
        Kcat_function(x, clearance_rate_constant = stats::coefficients(fitKcat))
      },
      colour = "blue"
    ) +
    ggplot2::annotate(
      "text",
      y = max(clear_curve_xy$y) * 0.8,
      x = max(clear_curve_xy$x) * 0.8,
      label = paste("R\u00b2=", r2),
      size = 5
    )

  print(diag_plot)

  #Make table with fit Kcat
  fit_95conf <- stats::confint(fitKcat)
  result <- data.frame(
    parameter = "rate_constant",
    estimate = stats::coefficients(fitKcat)[["clearance_rate_constant"]],
    lower = fit_95conf[[1]],
    upper = fit_95conf[[2]]
  )

  if (r2 < 0.8) {
    cli::cli_warn(
      "Poor fit of the depletion curve (R-squared {r2}): the data may not
       follow a one-phase decay."
    )
  } else if (min(clear_curve_xy$y) > 0.2 * max(clear_curve_xy$y)) {
    cli::cli_warn(
      "Little depletion: the rate constant cannot be estimated accurately."
    )
  }

  if (verbose) {
    .print_ivive_result(
      "fit_depletion_curve",
      inputs = list(data = data),
      result = result
    )
  }

  result
}
