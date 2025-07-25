% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm_fu_invitro.R
\name{Turner}
\alias{Turner}
\alias{calculate_fu_turner}
\title{Turner algorithm for Fu calculation}
\usage{
calculate_fu_turner(
  ionization,
  pka,
  ion_factor,
  log_lipophilicity,
  concentration_microsomes
)
}
\arguments{
\item{ionization}{Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)}

\item{pka}{vector of length of 2 with pkA of the compound}

\item{ion_factor}{ionization factor}

\item{log_lipophilicity}{LogP or LogMA of the compound}

\item{concentration_microsomes}{concentration of microsomes (in mg/mL)}

\item{blood_plasma_ratio}{Blood plasma ratio, this parameter is needed for Rodgers and Rowland and Poulin method for basic chemicals}

\item{fraction_unbound}{In Vivo Fraction Unbound in plasma from literature}
}
\value{
fuInvitro

fuInvitro
}
\description{
The different algorithms from literature to calculate Fu_mic

The Turner algorithm for calculating Fu in vitro
}
\details{
Algorithm4Fumic
}
\examples{
Turner("acid",3,100,3,1)
Halifax("base",3,100,3,1)
Halifax(c("base",0),c(3,0),100,3,1)
Austin("base",3,100,3,1)
Poulin("base",pka=3,X=100,Y=120,blood_plasma_ratio=1,fraction_unbound=0.2,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3,concentration_microsomes=1)
Poulin("acid",pka=3,X=100,concentration_cell_neutral_lipids=0.03,log_lipophilicity=3,concentration_microsomes=1)

}
