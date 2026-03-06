#' IVIVE for clearance
#' 
#' @description
#' IVIVE for clearance based on different type of values
#'
#' @param typeValue describe what type of in vitro data for clearance it is
#' @param units this are the units of the value. it does not matter typeValue is "decay experimentalcurve"
#' @param typeSystem microsomes or hepatocytes
#' @param expData depending on the type of value this can be a csv file or just one value
#' @param fu_invitro Fraction unbound in the in vitro hepatic system. if not known code will calculate it
#' @param empirical_scalar this is an option to include an extra empirical correction factor. Currently we are considering the scale factor of Wood 2017
#' @param tissue tissue of interest, since scaling factors are dependent on the tissue, will default to liver
#' @param species values can human and rat for now, defaulting to human
#' @param volMedium_mL volume of medium in the well (in mL)
#' @param REF relative expression or activity factor
#' @param cMicro_mgml concentration of microsome protein (in mg/mL)
#' @param cCells_Mml concentration of hepatocytes used (in million cells/mL)
#'
#' @return Specific clearance parameter (/min) to plug in PK-Sim
#' @export
#' @examples
#' exp_path<-system.file("extdata","clearance.csv",package="esqIVIVE")
#' expData<-read.csv(exp_path)
#' expData<-read.csv("tests/testthat/data/clearance.csv",header=TRUE)
#' IVIVE_clearance(typeValue="kcat",typeSystem="hepatocytes",cCells_Mml=0.5,units="/minutes",
#'                expData=2,fu_invitro=0.5,empirical_scalar="No")
#'
#' # example hepatocytes
#' IVIVE_clearance(typeValue="in vitro clearance parameter",typeSystem="hepatocytes",species="human",units="mL/minutes/Millioncells",
#' expData=18.27,fu_invitro=0.5,cCells_Mml=0.5,empirical_scalar="No")
#'
#' # example microsomes
#' IVIVE_clearance(typeValue="in vitro clearance parameter",typeSystem="microsomes",units="L/minutes/mg protein",
#'                 expData=18.27,fu_invitro=0.5,cMicro_mgml=0.5,volMedium_mL=0.5,empirical_scalar="No")
#'

IVIVE_clearance <- function(
  typeValue,
  units,
  expData,
  typeSystem,
  fu_invitro=1 ,
  empirical_scalar="No",
  tissue = "liver",
  species = "human",
  volMedium_mL = 1,
  REF = 1,
  cMicro_mgml = NULL,
  cCells_Mml = NULL
) {
  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes", "microsomes"))

  rlang::arg_match(
    typeValue,
    c("kcat", "halfLife", "in vitro clearance parameter")
  )

  #make it accordingly to microplate
  if (is.null(volMedium_mL)) {
    volMedium_mL <- 1
  } else {}

  if (is.null(tissue)) {
    tissue <- "liver"
  } else {}

  #Scaling factors

  scaling_factors <- list()
  scaling_factors[["human"]] <- data.frame(
    "organs" = c("brain", "gonads", "heart", "kidney", "gut", "liver", "lung"),
    "fcell" = c(0.96, 0.88, 0.76, 0.57, 0.88, 0.67, 0.23),
    "weightorgankgBW" = c(22.5, 0.167, 5.5, 6.67, 18.5, 32, 15),
    "MPGO" = rep(36, 7), #mg microsomal protein/g liver
    "HGO" = rep(119, 7)
  ) #million cells/g liver

  scaling_factors[["rat"]] <- data.frame(
    "organs" = c("brain", "gonads", "heart", "kidney", "gut", "liver", "lung"),
    "fcell" = c(0.96, 0.79, 0.64, 0.7, 0.88, 0.72, 0.19),
    "weightorgankgBW" = c(7.39, 10, 3.47, 10, 31, 43, 4.34),
    "MPGO" = rep(50, 7), #mg microsomal protein/g liver
    "HGO" = rep(122, 7)
  ) #million cells/g liver

  scaling_factors[["dog"]] <- data.frame(
    "organs" = c("brain", "gonads", "heart", "kidney", "gut", "liver", "lung"),
    "fcell" = c(0.96, 0.79, 0.64, 0.7, 0.88, 0.72, 0.19), #factor of the rat
    "weightorgankgBW" = c(5.88, 0.36, 5.1, 11.7, 17.5, 25, 10.3),
    "MPGO" = rep(50, 7), #mg microsomal protein/g liver
    "HGO" = rep(201, 7)
  ) #million cells/g liver

  #Get species scaling factor
  SF_matrix <- scaling_factors[[species]]
  fintcell <- SF_matrix[which(SF_matrix[, 1] == tissue), "fcell"]
  liverkgBW <- SF_matrix[which(SF_matrix[, 1] == tissue), "weightorgankgBW"] #is only used in some cases

  #chose the system specific scaling factors
  if (typeSystem == "microsomes") {
    nLiver <- SF_matrix[which(SF_matrix[, 1] == tissue), "MPGO"] # mg/g liver
    cInvitro <- cMicro_mgml #mg/mL
  } else if (typeSystem == "hepatocytes") {
    nLiver <- SF_matrix[which(SF_matrix[, 1] == tissue), "HGO"]
    cInvitro <- cCells_Mml # million cells/mL assay
  } else {}

  #combination scale factors
  SF <- nLiver / fintcell * REF / fu_invitro
  SF2 <- 1 / fintcell * REF / fu_invitro

  #Derive the in vitro clearance value---------------------
  if (typeValue == "kcat") {
    
    kcat_min <- expData
    ClspePermin <- kcat_min / cInvitro * SF
    
  } else if (typeValue == "halfLife") {
    halfLife <- expData

    #make matrix of calculations depedning on the units
    matrixCalcHalf <- cbind(c("minutes", "hours", "seconds"), c(1, 1 / 60, 60))

    multFactorHalf <- as.double(matrixCalcHalf[
      which(matrixCalcHalf[, 1] == units),
      2
    ])
    kcat_min <- 0.693 / halfLife * multFactorHalf
    ClspePermin <- kcat_min / cInvitro * SF
  } else if (typeValue == "in vitro clearance parameter") {
    #make matrix of calculations depending on the units
    matrixCalClear <- cbind(
      c(
        "/minutes",
        "/hours",
        "/seconds",
        "mL/minutes/millioncells",
        "uL/minutes/millioncells",
        "L/minutes/millioncells",
        "mL/hours/millioncells",
        "uL/hours/millioncells",
        "L/hours/millioncells",
        "mL/seconds/millioncells",
        "uL/seconds/millioncells",
        "L/seconds/millioncells",
        "mL/minutes/cell",
        "uL/minutes/cell",
        "mL/hours/cell",
        "uL/hours/cell",
        "mL/seconds/cell",
        "uL/seconds/cell",
        "mL/minutes/mg protein",
        "uL/minutes/mg protein",
        "L/minutes/mg protein",
        "mL/hours/mg protein",
        "uL/hours/mg protein",
        "L/hours/mg protein",
        "mL/seconds/mg protein",
        "uL/seconds/mg protein",
        "L/seconds/mg protein",
        "mL/minutes",
        "uL/minutes",
        "mL/seconds",
        "uL/seconds",
        "mL/hours",
        "uL/hours",
        "mL/minutes/kg",
        "uL/minutes/kg",
        "mL/hours/kg",
        "uL/hours/kg"
      ),
      c(
        1 / cInvitro,
        1 / cInvitro * 60,
        1 / cInvitro / 60,
        1,
        1 / 1000,
        1000,
        1 / 60,
        1 / 60000,
        1000 / 60,
        60,
        60 / 1000,
        60000,
        1E6,
        1000,
        1E6 / 60,
        1E6 / 60000,
        60 * 1E6,
        60 * 1E6 / 1000,
        1,
        1 / 1000,
        1000,
        1 / 60,
        1 / 60000,
        1000 / 60,
        60,
        60 / 1000,
        60000,
        1 / (cInvitro * volMedium_mL),
        1 / (cInvitro * volMedium_mL) / 1000,
        1 / (cInvitro * volMedium_mL) * 60,
        1 / (cInvitro * volMedium_mL) * 60 / 1000,
        1 / (cInvitro * volMedium_mL) / 60,
        1 / (cInvitro * volMedium_mL) / 60000,
        1 / liverkgBW,
        1 / liverkgBW / 1000,
        1 / liverkgBW / 60,
        1 / liverkgBW / 60000
      )
    ) #multiplying by the g liver weight per kilo BW

    matrixCalClear[1:33, 2] <- as.double(matrixCalClear[1:33, 2]) * SF
    matrixCalClear[34:37, 2] <- as.double(matrixCalClear[34:37, 2]) * SF2

    multFactorCl <- as.double(matrixCalClear[
      which(matrixCalClear[, 1] == units),
      2
    ])

    if (length(multFactorCl) == 0) {
      print("Not identified units")
    } else {}

    ClspePermin <- expData * multFactorCl
  } else {
    ClspePermin <- 0
  }

  #Add scalars from Wood et al 2017-https://doi.org/10.1124/dmd.117.077040.
  wood_sf <- list()
  wood_sf[["human"]] <- data.frame(
    Cl_ranges = c("<10", "10-100", "100-1000", "1000-1000", ">10000"),
    microsomes = c(0.7, 1.8, 4.6, 7.5, 58),
    hepatocytes = c(0.61, 3.9, 7.1, 22, 1200)
  )
  wood_sf[["rat"]] <- data.frame(
    Cl_ranges = c("<10", "10-100", "100-1000", "1000-1000", ">10000"),
    microsomes = c(0.086, 0.83, 1.7, 2.5, 230),
    hepatocytes = c(0.13, 1.6, 3.2, 7.2, 180)
  )
  wood_table <- wood_sf[[species]]

  if (empirical_scalar == "Yes") {
    # conditional to pick the right wood scalar
    # this scalar are empitical and they are based on in vitro data tending to overestimate slow clearance
    # and underpredict fast clearance
    if (ClspePermin * liverkgBW < 10) {
      ClspePermin <- ClspePermin *
        wood_table[
          which(wood_table[, 1] == "<10"),
          which(colnames(wood_table) == typeSystem)
        ]
    } else if (ClspePermin * liverkgBW < 100 && ClspePermin * liverkgBW > 10) {
      ClspePermin <- ClspePermin *
        wood_table[
          which(wood_table[, 1] == "10-100"),
          which(colnames(wood_table) == typeSystem)
        ]
    } else if (
      ClspePermin * liverkgBW < 1000 && ClspePermin * liverkgBW > 100
    ) {
      ClspePermin <- ClspePermin *
        wood_table[
          which(wood_table[, 1] == "100-1000"),
          which(colnames(wood_table) == typeSystem)
        ]
    } else if (
      ClspePermin * liverkgBW < 10000 && ClspePermin * liverkgBW > 1000
    ) {
      ClspePermin <- ClspePermin *
        wood_table[
          which(wood_table[, 1] == "1000-10000"),
          which(colnames(wood_table) == typeSystem)
        ]
    } else if (ClspePermin * liverkgBW > 10000) {
      ClspePermin <- ClspePermin *
        wood_table[
          which(wood_table[, 1] == ">10000"),
          which(colnames(wood_table) == typeSystem)
        ]
    }
  } else {}

  return("ClspePermin" = as.double(ClspePermin))
}
