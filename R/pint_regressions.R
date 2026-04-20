#' Function to convert from Papp to Pint transcellular permeability from a empirical regression
#' 
#' @param Papp_cms , permeability from Caco-2 in cms
#'
#' @returns pint in units cm/s
#' @export
#'
#' @examples
#'pint_caco2_empir(Papp_cms=2.3E-6)
#'  
pint_caco2_empir<-function(Papp_cms){
  #based on calibration with high solubility
  pint_cms<-0.0001*10^((0.4428*(log10(Papp_cms*1000000)-3.0941)))
  return(pint_cms)
  
}


#' Function to convert from Papp to Pint transcellular permeability from a empirical regression
#'
#' @param Peff_cms Peff as obtained in vivo or QSARs in cms
#'
#' @returns pint in units cm/s
#' @export
#'
#' @examples
#' pint_peff_empir(Peff_cms=2.3E-6)
#' 
pint_peff_empir<-function(Peff_cms){
  #based on calibration with high solubility
  pint_cms<-0.0001*10^((0.864*(log10(Peff_cms*1000)-3.1029)))
  return(pint_cms)
  
}