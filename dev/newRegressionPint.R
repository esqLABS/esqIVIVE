#' Function to convert from Papp to Pint transcellular permeability from a empirical regression
#' 
#' @param Papp_cms , permeability from Caco-2 in cms
#'
#' @returns pint in units cm/s
#' @export
#'
#' @examples
#' 
#feed table with permeability (specific units, for the specifc chemicals)
newRegressionPint<-function(table_expperm_cms,ionization="considered"){
  #pull fitted Pint
  #bind tables
  #identify if there is enough coverage, 3 chemical per low, moderate and high perm
  #experimental
  #linear fit for the log
  
  
  
  return(list("slope"=slope,"y0"=y0,"R2"=R2))
}