\name{scenariosFromStatistics}
\alias{scenariosFromStatistics}

\title{Make Scenarios to Match Capture Statistics}

\description{

  The \code{\link{make.scenarios}} function requires prior knowledge of
  population density and the intercept of the detection function
  (g0). This function provides an alternative mechanism for generating
  scenarios from a value of sigma and target values for the numbers of
  individuals n and recaptures r. Only a halfnormal detection function
  is supported (probability, not hazard), and many options in
  \code{\link{make.scenarios}} have yet to be implemented. Only a single
  detector layout and single mask may be specified.
  
}

\usage{

scenariosFromStatistics(sigma, noccasions, traps, mask, nval, rval,
 g0.int = c(0.001, 0.999))

}

\arguments{
  \item{sigma}{numeric vector of one or more values for sigma}
  \item{noccasions}{integer vector of number of sampling occasions}
  \item{traps}{traps object}
  \item{mask}{mask object}
  \item{nval}{integer vector of values of n}
  \item{rval}{integer vector of values of r}
  \item{g0.int}{numeric vector defining the interval to be searched for g0 }
}
\details{

  The algorithm is based on R code in Appendix B of Efford, Dawson and
  Borchers (2009).

}

\value{

  A scenario dataframe with one row for each combination of \code{sigma},
  \code{noccasions}, \code{nval} and \code{rval}.
  
}

\references{
  
Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\seealso{

  \code{\link{make.scenarios}}
  
}

\examples{

grid36 <- make.grid(nx = 6, ny = 6, spacing = 200)
mask <- make.mask(grid36, buffer = 2000)
scen <- scenariosFromStatistics (sigma = c(200,400), noccasions = 44,
    traps = grid36, mask = mask, nval = 14, rval = 34)
sim <- run.scenarios(scen, nrepl = 5, traps = grid36, mask = mask)
summary(sim)

}

\keyword{ datagen }
