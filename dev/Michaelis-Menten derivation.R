#Code to derive convert Michaelis-Menten
#' Title
#'
#' @param typeValue typeValue, default is Vmax&Km but can be exp_curve
#' @param typeSystem if hepatocytes or microsomes
#' @param fu_invitro value of f 
#' @param vmax 
#' @param km_micromolar 
#' @param expData 
#' @param tissue 
#' @param species 
#' @param volume_medium 
#' @param REF 
#' @param concentration_microsomes 
#' @param concentration_cells 
#'
#' @returns Vmax and Km
#' @export
#'
#' @examples


IVIVE_MM <- function(
  typeValue=Vmax&Km,
  typeSystem,
  fu_invitro,
  vmax = NULL,
  km_micromolar = NULL,
  expData = NULL,
  tissue = NULL,
  species = NULL,
  volume_medium = NULL,
  REF = NULL,
  concentration_microsomes = NULL,
  concentration_cells = NULL
) {
  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes", "microsomes"))

  # check if the arguments are valid
  rlang::arg_match(typeValue, c("exp_curve", "Vmax&Km"))

  if (is.null(fu_invitro)) {
    fu_invitro <- 1
  } else if (fu_invitro == 0) {
    print("problem fu_invitro=0")
  } else {}
  if (is.null(species)) {
    species <- "human"
  } else {}

  #make it accordingly to microplate
  if (is.null(volume_medium)) {
    volume_medium <- 1
  } else {}

  if (is.null(tissue)) {
    tissue <- "liver"
  } else {}

  if ("typeValue" == exp_curve) {
    MM <- get_MM("mm_exp_curve_xy" = mm_exp_curve)
    km_micromolar <- MM[1, 1]
    vmax <- MM[2, 1]
    
  } else {
    
    km_micromolar <- km_micromolar
    vmax <- vmax
  }

  #Calculate fraction unbound ---------------------------------------------------

  #Correct Km for fraction unbound

  Km_unb <- Km * fu_invitro

  #Calculate in vivo Cl_u--------------------------------------------------------
  sf_cyp = 108 #pmol_permgmicro, specific of enzyme

  #Calculate the specific clearance, it will be in units of per min
  if (invitro_system == "microsomes") {
    N_liver = N_liver = 40 # mg/g
  } else if (invitro_system == "hepatocytes") {
    N_liver = 139 * 10**6 # cell/g
  } else {
    print("Error in vitro system")
  }

  fu_intracell_liver = 0.67

  return(list(Vmax, Km_unb_uM))
}

