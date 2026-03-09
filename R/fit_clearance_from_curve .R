
#' function to determine clearance from experimental curve
#' @name fit_clearance_from_curve
#' @param expData_tmin_cuM table with first column as time in minutes and second column as concetration in uM
#' 
#'
#' @returns kcat in per min, this is not clearance ready for pksim
#' @export
#'
#' @examples
#' exp_path<-system.file("extdata","clearance.csv",package="esqIVIVE")
#' expData<-read.csv(exp_path)
#' fit_clearance_from_curve(expData)
#' 
#' 
fit_clearance_from_curve <- function(expData_tmin_cuM) {
  library(ggplot2)
  #Load the depletion curve
  clear_curve_xy <- expData_tmin_cuM
  colnames(clear_curve_xy) <- c("x", "y")
  
  #find the starting concentration
  y0 = mean(clear_curve_xy$y[clear_curve_xy$x == 0])
  
  #create clearance model
  Kcat_function <- function(x, clearance_rate_constant) {
    y = y0 * exp(-clearance_rate_constant * x)
    return(y)
  }
  
  #fit model
  fitKcat <- nls(
    y ~ Kcat_function(x, clearance_rate_constant),
    data = clear_curve_xy,
    start = list(clearance_rate_constant = 0.01),
    trace = TRUE
  )
  
  r_squared_nls <- function(model) {
    rss <- sum(residuals(model)^2)
    tss <- sum((clear_curve_xy$x- 
                  mean(clear_curve_xy$x))^2)  # y - mean(y)
    1 - (rss / tss)
  }
  r2<-round(r_squared_nls(fitKcat),digits=3)
  
  colnames(clear_curve_xy) <- colnames(expData)
  #Plot for evaluating if fit is reasonable
  diag_plot <- ggplot(clear_curve_xy, aes(x = Time_min, y = Concentration_uM)) +
    geom_point() +
    labs(title = "fit curve") +
    stat_function(
      fun = function(x) Kcat_function(x, clearance_rate_constant = coefficients(fitKcat)),
      colour = "blue")+
    annotate("text", y = max(clear_curve_xy$Concentration_uM)*0.8,
             x = max(clear_curve_xy$Time_min)*0.8, 
             label = paste("R²=",r2), size = 5)
  
  print(diag_plot)
  
  #Make table with fit Kcat
  fit_95conf = confint(fitKcat)
  kcat = as.double(data.frame(
    "Kcat_permin" = coefficients(fitKcat)["clearance_rate_constant"],
    fit_95conf[1],
    fit_95conf[2]
  ))
  names(kcat) = c("Mean_min-1", "2.5%_CI", "95%_CI")
  
  if (r2<0.8){
    warning("poor fit of clearance curve")
  } else if (min(clear_curve_xy$Concentration_uM)>0.2*max(clear_curve_xy$Concentration_uM)){
    warning("little depletion, cannot be fit accuratly ")
  } else {}
  
  return(kcat)
}
