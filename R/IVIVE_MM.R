#Code to derive convert Michaelis-Menten
#' Title
#' @name IVIVE_MM
#' @description function that scales Vmax and correct km for fraction unbound
#' 
#' @param typeSystem if hepatocytes or microsomes
#' @param fu_invitro value of fractionunbound in vitro, the default is 1
#' @param vmax as umol/min/million hepatocytes or umol/min/mg microsomal protein
#' @param km_micromolar Km of the enzyme reaction, in uM
#' @param tissue liver, brain, lung, kidney, gonads and gut, default is liver
#' @param species human, rat or dog, default human
#' @param REF relative expression or activity factor, default is 1.
#' To use this option the reference concentration of
#'  of the enzyme of interest in pksim needs to be 1 uM
#' @param verbose if TRUE, print the inputs and the resulting Vmax/Km
#' @returns Vmax in umol/min/L and Km_unb in uM
#' @export
#'
#' @examples
#' IVIVE_MM (typeSystem="hepatocytes",vmax=2,km_micromolar=1,tissue="liver",species="human",REF=1)
#' IVIVE_MM (typeSystem="microsomes",fu_invitro=0.2,vmax=2,km_micromolar=1)

IVIVE_MM <- function(
  typeSystem,
  fu_invitro=1,
  vmax,
  km_micromolar,
  tissue =  "liver",
  species = "human",
  REF = 1,
  verbose = FALSE
) {
  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes", "microsomes"))

  if (fu_invitro == 0) {
    print("problem fu_invitro=0")
  } else if (fu_invitro >1 ){
    print("problem fu_invitro>1")
  }

  #Correct Km for fraction unbound
  Km_unb_uM <- km_micromolar * fu_invitro

  #Calculate in vivo Vmax--------------------------------------------------------
  #Scaling factors
  path <- system.file("extdata", "scaling_factors.csv", package = "esqIVIVE")
  scaling_factors<-read.csv(path)

  # check if the arguments are valid
  rlang::arg_match(species, unique(scaling_factors[,"species"]))
  rlang::arg_match(tissue, unique(scaling_factors[,"organ"]))

  #Get species scaling factor
  species_row<-which(scaling_factors[,"species"]==species)
  organ_row<-which(scaling_factors[,"organ"]==tissue)
  overlap_row<-intersect(species_row,organ_row)
  fintcell <- scaling_factors[overlap_row, "fcell"]
 
  #chose the system specific scaling factors
  if (typeSystem == "microsomes") {
    scfactor <- scaling_factors[overlap_row, "MicProtGO"] # mg protein/g liver
 } else if (typeSystem == "hepatocytes") {
   scfactor <- scaling_factors[overlap_row, "CellsGO"]
 } else {
    warning("typeSystem not identified, only hepatocytes or microsomes allowed")
 }
  dens<-1000 #g/L
  vmax_umol_minL  <- vmax*scfactor*REF/fintcell*dens

  result <- list("vmax_umol_minL"=vmax_umol_minL,"Km_unb_uM"=Km_unb_uM)

  if (verbose) {
    .print_ivive_result(
      "IVIVE_MM",
      inputs = list(typeSystem = typeSystem, fu_invitro = fu_invitro, vmax = vmax, km_micromolar = km_micromolar, tissue = tissue, species = species, REF = REF),
      result = result
    )
  }

  return(result)
}

