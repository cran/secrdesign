\name{costing}
\alias{costing}

\title{Cost of SECR design}

\description{

The cost of implementing a spatially explicit capture--recapture design depends on the detector layout, the number of detections and the various unit costs.

}

\usage{
costing(traps, nr, noccasions, unitcost = list(), nrepeats = 1, routelength = NULL, 
    setupoccasion = TRUE)
}

\arguments{
  \item{traps}{traps object for detector array}
  \item{nr}{numeric vector with E(\eqn{n}) and E(\eqn{r}) as first two elements}
  \item{noccasions}{integer number of sampling occasions}
  \item{unitcost}{list with unit costs (see Details)}
  \item{nrepeats}{integer number of repeated arrays}
  \item{routelength}{numeric route length (km)}
  \item{setupoccasion}{logical; if TRUE then the cost of a setup visit is included (noccasions+1) }
}


\details{

\code{nr} is a vector with the expected sample sizes (numbers of individuals and recaptures), usually the output from \code{\link{Enrm}}.

\code{unitcost} should be a list with at least one of the components `perkm', `perarray', `perdetector', `pervisit' and `perdetection'. 

The number of occasions (noccasions) is incremented by 1 if \code{setupoccasion} is TRUE.

\tabular{llll}{
Component \tab Unit cost \tab Costing \cr
Arrays \tab perarray \tab perarray x nrepeats \cr
Detectors \tab perdetector \tab perdetector x nrow(traps) x nrepeats \cr
Travel \tab  perkm \tab perkm x routelength x noccasions x nrepeats \cr
Visits \tab pervisit \tab sum(pervisit x trapcost) x noccasions x nrepeats \cr
Detections \tab perdetection \tab perdetection x total detections (E(\eqn{n}) + E(\eqn{r})) \cr
}

`Travel' and `Visits' are alternative ways to cost field time. The variable `routelength' represents the length of a path followed to visit all detectors; if not specified it is approximated by the sum of the nearest-trap distances. The variable `trapcost' is a vector of length equal to the number of detectors. By default it is a vector of 1's, but detector- specific values may be provided as trap covariate `costpervisit'. In the latter case the value of `pervisit' should probably be 1.0.

`Arrays' and `Detectors` represent one-off costs.

`Detections' includes costs such as handling time and laboratory DNA analysis.

See  \url{../doc/secrdesign-tools.pdf} for more.

}

\value{

A named numeric vector

}

\seealso{

\code{\link{Enrm}}, \code{\link{scenarioSummary}}

}

\examples{

tr <- make.grid(8, 8, spacing = 25)
msk <- make.mask(tr, buffer = 100, type = 'trapbuffer')
nrm <- Enrm(D = 5, tr, msk, list(lambda0 = 0.2, sigma = 20), 5)
costing (tr, nrm, 5, unitcost = list(pervisit = 5, perdetection = 15))

}
