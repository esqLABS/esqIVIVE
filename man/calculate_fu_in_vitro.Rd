% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_fu_inc_mic.R
\name{calculate_fu_in_vitro}
\alias{calculate_fu_in_vitro}
\title{getInVitroFractionUnbound}
\usage{
calculate_fu_in_vitro(
  partition_qspr,
  log_lipophilicity,
  ionization,
  type_system,
  fetal_bovine_serum_fraction,
  microplate_type,
  volume_medium,
  pka = NULL,
  henry_law_constant = NULL,
  fraction_unbound = NULL,
  blood_plasma_ratio = NULL,
  concentration_microsomes = NULL,
  concentration_cells = NULL
)
}
\arguments{
\item{partition_qspr}{type of assumption used (Poulin and Theil, PK-Sim® Standard, Rodgers & Rowland, Schmidtt,
then from literature, Poulin, Turner, Austin and Halifax. See thevignette for more details)}

\item{log_lipophilicity}{LogP or LogMA of the compound}

\item{ionization}{Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)}

\item{type_system}{microsomes or hepatocytes}

\item{fetal_bovine_serum_fraction}{fraction of serum concentration, values can only go from 0-1}

\item{microplate_type}{number of wells in the microplate}

\item{volume_medium}{volume of medium in the well (in mL)}

\item{pka}{vector of length of 2 with pkA of the compound}

\item{henry_law_constant}{Henry's Law Constant (in atm/(m3*mol))}

\item{fraction_unbound}{In Vivo Fraction Unbound in plasma from literature}

\item{blood_plasma_ratio}{Blood plasma ratio, this parameter is needed for Rodgers and Rowland and Poulin method for basic chemicals}

\item{concentration_microsomes}{concentration of microsomes (in mg/mL)}

\item{concentration_cells}{concentration of cells (in million cells/mL)}
}
\value{
fuInvitro and possible warning for evaporation
}
\description{
Compute the fraction unbound in vitro
}
\details{
mayeb consider to have average data..
}
\examples{
calculate_fu_in_vitro(
 partition_qspr = "All PK-Sim Standard", log_lipophilicity = 3, ionization = c("acid", 0),
 type_system = "hepatocytes", fetal_bovine_serum_fraction = 0, microplate_type = 96,
 volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2)

calculate_fu_in_vitro(
 partition_qspr = "Poulin and Theil + fu", log_lipophilicity = 3, ionization = c("acid", 0),
 type_system = "hepatocytes", fetal_bovine_serum_fraction = 0, microplate_type = 96,fraction_unbound=0.01,blood_plasma_ratio=2,
 volume_medium = 0.22, pka = c(6, 0), henry_law_constant = 1E-6, concentration_cells = 2)

calculate_fu_in_vitro(
 partition_qspr = "All Schmitt", log_lipophilicity = 0.42, ionization = c("acid", 0),
 type_system = "microsomes", fetal_bovine_serum_fraction = 0, microplate_type = 96,fraction_unbound=0.2,blood_plasma_ratio=1,
 volume_medium = 0.22, pka = c(6, 0), concentration_microsomes = 2)

}
