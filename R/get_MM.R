#' Get michaelis-menten  parameters form experimental curves
#' @name get_MM
#'
#' @description
#' Function to derive micahelis-menten form raw data
#'
#' @param experimental_conc_velocity is a experimental curve with concentration in the first column and velocity in the second column
#'
#' @return fitresults_vmax_km
#' @examples
#' mm_curve_path<-system.file("extdata","michaelis_menten_curve.csv","esqIVIVE")
#' mm_curve<-read.csv(mm_curve_path)
#' get_MM(mm_curve)
#' 
#' @export
#' 
#' #function to get MM form raw data
  get_MM <- function(experimental_conc_velocity) {
    library(ggplot2)
    colnames(experimental_conc_velocity)<-c("Concentration","Velocity")
    
    #fit model
    fitmm <- nls(
      Velocity ~ Vmax * Concentration / (Km + Concentration),
      data = experimental_conc_velocity,
      start = list(
        Vmax = max(experimental_conc_velocity$Velocity),
        Km = mean(experimental_conc_velocity$Concentration)
      ),
      trace = FALSE
    )
    fit_95conf = confint(fitmm)
    
    r_squared_nls <- function(model) {
      rss <- sum(residuals(model)^2)
      tss <- sum((experimental_conc_velocity$Velocity- 
                  mean(experimental_conc_velocity$Velocity))^2)  # y - mean(y)
      1 - (rss / tss)
    }
    r2<-round(r_squared_nls(fitmm),digits=3)
    
    #check if fitting is good
    mm_fuction <- function(Concentration) {
      Velocity = coefficients(fitmm)["Vmax"] * Concentration / 
                (coefficients(fitmm)["Km"] + Concentration)
      return(Velocity)
    }
    
    plot_diagnosis <- ggplot(
      experimental_conc_velocity,
      aes(x = Concentration, y = Velocity)
    ) +
      geom_point() +
      stat_function(fun = function(x) mm_fuction(x), colour = "blue")+
      annotate("text", x = max(experimental_conc_velocity$Concentration)*0.8,
               y = max(experimental_conc_velocity$Velocity)*0.8, 
               label = paste("R²=",r2), size = 5)
    
    #Add fit values in dataframe for calculations
    fitresults_vmax_km = data.frame(
      "Mean" = c(coefficients(fitmm)["Km"], coefficients(fitmm)["Vmax"]),
      "2.5_percent" = c(fit_95conf["Km", 1], fit_95conf["Vmax", 1]),
      "95%_percent" = c(fit_95conf["Km", 2], fit_95conf["Vmax", 2])
    )
    
    row.names(fitresults_vmax_km) = c("Km_uM", "Vmax_umol_min_mgmicroORcells")
    print(plot_diagnosis)
    return(fitresults_vmax_km)
 }
