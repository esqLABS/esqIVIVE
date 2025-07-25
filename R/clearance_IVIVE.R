#' IVIVE for clearance
#'
#' @description
#' IVIVE for clearance
#'
#' @param typeValue describe what type of in vitro data for clearance it is
#' @param units this are the units of the value. it does not matter typeValue is "decay experimentalcurve"
#' @param typeSystem microsomes or hepatocytes
#' @param expData depending on the type of value this can be a csv file or just one value
#' @param fu_invitro Fraction unbound in the in vitro hepatic system. if not known code will calculate it
#' @param empirical_scalar this is an option to include an extra empirical correction factor. Currently we are considering the scale factor of Wood 2017
#' @param tissue tissue of interest, since scaling factors are dependent on the tissue, will default to liver
#' @param species values can human and rat for now, defaulting to human
#' @param volMedium volume of medium in the well (in mL)
#' @param REF relative expression or activity factor
#' @param cMicro_mgml concentration of microsome protein (in mg/mL)
#' @param cCells_Mml concentration of hepatocytes used (in million cells/mL)
#' @param partition_qspr type of assumption used for calculating fu_invitro (optional)
#' @param log_lipophilicity LogP or LogMA of the compound (optional)
#' @param ionization vector of length 2 with ionization class, acid, neutral and base (optional)
#' @param fetal_bovine_serum_fraction fraction of serum concentration, values can only go from 0-1 (optional)
#' @param microplate_type number of wells in the microplate (optional)
#' @param volume_medium volume of medium in the well (in mL) (optional)
#' @param pka vector of length 2 with pKa values of the compound (optional)
#' @param henry_law_constant Henry's Law Constant (in atm/(m3*mol)) (optional)
#' @param fraction_unbound In Vivo Fraction Unbound in plasma from literature (optional)
#' @param blood_plasma_ratio Blood plasma ratio (optional)
#'
#' @return Specific clearance parameter to plug in PK-Sim
#' @export
#' @examples
#'
#' expData<-read.csv("tests/testthat/data/clearance.csv",header=TRUE)
#' clearance_IVIVE(typeValue="decay experimentalcurve",typeSystem="hepatocytes",cCells_Mml=0.5,units="/minutes",
#'                expData=expData,fu_invitro=0.5,empirical_scalar="No")
#'
#' # example hepatocytes
#' clearance_IVIVE(typeValue="in vitro clearance parameter",typeSystem="hepatocytes",species="human",units="mL/minutes/Millioncells",
#' expData=18.27,fu_invitro=0.5,cCells_Mml=0.5,empirical_scalar="No")
#'
#' # example microsomes
#' predict_clearance_ivive(typeValue="in vitro clearance parameter",typeSystem="microsomes",units="L/minutes/mg protein",
#'                 expData=18.27,fu_invitro=0.5,cMicro_mgml=0.5,volMedium=0.5,empirical_scalar="No")
#'

predict_clearance_ivive <- function(
  typeValue,
  units,
  expData,
  typeSystem,
  fu_invitro = NULL,
  empirical_scalar = "No",
  tissue = NULL,
  species = NULL,
  volMedium = NULL,
  REF = NULL,
  cMicro_mgml = NULL,
  cCells_Mml = NULL,
  # Additional parameters for calculating fu_invitro
  partition_qspr = NULL,
  log_lipophilicity = NULL,
  ionization = NULL,
  fetal_bovine_serum_fraction = NULL,
  microplate_type = NULL,
  volume_medium = NULL,
  pka = NULL,
  henry_law_constant = NULL,
  fraction_unbound = NULL,
  blood_plasma_ratio = NULL
) {
  # check if the arguments are valid
  rlang::arg_match(typeSystem, c("hepatocytes", "microsomes"))

  rlang::arg_match(
    typeValue,
    c("decay experimentalcurve", "halfLife", "in vitro clearance parameter")
  )

  # Calculate fu_invitro if parameters are provided
  if (is.null(fu_invitro) && !is.null(partition_qspr)) {
    # Set default values for ionization and pka if not provided
    if (is.null(ionization)) {
      ionization <- c("neutral", 0)
    }
    if (is.null(pka)) {
      pka <- c(0, 0)
    }
    if (is.null(henry_law_constant)) {
      henry_law_constant <- 0.000001
    }

    fu_invitro <- calculate_fu_in_vitro(
      partition_qspr = partition_qspr,
      log_lipophilicity = log_lipophilicity,
      ionization = ionization,
      type_system = typeSystem,
      fetal_bovine_serum_fraction = fetal_bovine_serum_fraction,
      microplate_type = microplate_type,
      volume_medium = volume_medium,
      pka = pka,
      henry_law_constant = henry_law_constant,
      fraction_unbound = fraction_unbound,
      blood_plasma_ratio = blood_plasma_ratio,
      concentration_microsomes = cMicro_mgml,
      concentration_cells = cCells_Mml
    )
  }

  #defaults for empty variable
  if (is.null(fu_invitro)) {
    fu_invitro <- 1
  } else if (fu_invitro == 0) {
    print("problem fu_invitro=0")
  } else {}

  if (is.null(REF)) {
    REF <- 1
  } else {}

  if (is.null(species)) {
    species <- "human"
  } else {}

  #make it accordingly to microplate
  if (is.null(volMedium)) {
    volMedium <- 1
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
  if (typeValue == "decay experimentalcurve") {
    kcat_min <- fit_clearance_from_curve(expData)[1]
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
        "mL/minutes/Millioncells",
        "uL/minutes/Millioncells",
        "L/minutes/Millioncells",
        "mL/hours/Millioncells",
        "uL/hours/Millioncells",
        "L/hours/Millioncells",
        "mL/seconds/Millioncells",
        "uL/seconds/Millioncells",
        "L/seconds/Millioncells",
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
        1 / (cInvitro * volMedium),
        1 / (cInvitro * volMedium) / 1000,
        1 / (cInvitro * volMedium) * 60,
        1 / (cInvitro * volMedium) * 60 / 1000,
        1 / (cInvitro * volMedium) / 60,
        1 / (cInvitro * volMedium) / 60000,
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

#function to determine clearance from experimental curve
fit_clearance_from_curve <- function(expData) {
  library(ggplot2)

  #Load the depletion curve
  clear_curve_xy <- expData
  colnames(clear_curve_xy) <- c("x", "y")

  #find the starting concentration
  y0 = mean(clear_curve_xy$y[clear_curve_xy$x == 0])

  #create clearance model
  Kcat_function <- function(initial_concentration, time, clearance_rate_constant) {
    y = initial_concentration * exp(-clearance_rate_constant * time)
    return(y)
  }

  #fit model
  fitKcat <- nls(
    y ~ Kcat_function(initial_concentration, x, clearance_rate_constant),
    data = clear_curve_xy,
    start = list(clearance_rate_constant = 0.01),
    trace = TRUE
  )

  colnames(clear_curve_xy) <- colnames(expData)
  #Plot for evaluating if fit is reasonable
  diag_plot <- ggplot(clear_curve_xy, aes(x = Time_min, y = Concentration_uM)) +
    geom_point() +
    labs(title = "fit curve") +
    stat_function(
      fun = function(x) Kcat_function(x, initial_concentration = y0, clearance_rate_constant = coefficients(fitKcat)),
      colour = "blue"
    )

  print(diag_plot)

  #Make table with fit Kcat
  fit_95conf = confint(fitKcat)
  kcat = as.double(data.frame(
    "Kcat_permin" = coefficients(fitKcat)["clearance_rate_constant"],
    fit_95conf[1],
    fit_95conf[2]
  ))
  names(kcat) = c("Mean", "2.5%_CI", "95%_CI")

  return(kcat)
}
