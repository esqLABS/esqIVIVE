rm(list = ls())

expData <- read.csv("tests/Clearance_example.csv", header = TRUE)
library(ggplot2)

source("R/clearance_IVIVE.R")

# example for curve
clearance_IVIVE(
  typeValue = "decay experimentalcurve", expData = expData,
  partitionQSPR = "All Schmitt",
  logLipo = 5,
  ionization = c("acid", 0, 0),
  typeSystem = "hepatocytes",
  FBS = 0.0, pKa = c(3, 0, 0),
  hlcAt = 1E-6,
  microplateType = 96,
  volMedium = 0.2,
  cCells = 0.02
)

# example for half life
clearance_IVIVE(
  typeValue = "halfLife", expData = 3,
  units = "hours",
  partitionQSPR = "All Schmitt",
  logLipo = 5,
  ionization = c("acid", 0, 0),
  typeSystem = "microsomes",
  FBS = 0, pKa = c(3, 0, 0),
  hlcAt = 1E-7,
  microplateType = 96,
  volMedium = 0.2,
  cMicro = 0.2
)

# example for clearance values
clearance_IVIVE(
  typeValue = "in vitro clearance parameter", expData = 0.3,
  units = "mL/seconds/mg protein",
  partitionQSPR = "All Poulin and Theil",
  logLipo = 2,
  ionization = c("neutral", 0, 0),
  typeSystem = "hepatocytes",
  FBS = 0, pKa = c(0, 0, 0),
  hlcAt = 0.02,
  microplateType = 96,
  volMedium = 0.2,
  cCells = 0.2
)


# test clearance with given fu_hep
clearance_IVIVE(
  typeSystem = "microsomes", typeValue = "in vitro clearance parameter", expData = 0.3,
  units = "uL/seconds/mg protein",
  partitionQSPR = "All Poulin and Theil",
  logLipo = 2,
  FBS = 0, fu_hep = 0.1,
  microplateType = 96,
  cMicro = 0.2
)
