\name{secrdesign-package}
\alias{secrdesign-package}
\alias{secrdesign}

\docType{package}

\title{Spatially Explicit Capture--Recapture Study Design }

\description{
Tools to assist the design of spatially explicit capture--recapture studies of animal populations.
}

\details{
\tabular{ll}{
Package: \tab secr\cr
Type: \tab Package\cr
Version: \tab 2.3.0\cr
Date: \tab 2015-03-25\cr
License: \tab GNU General Public License Version 2 or later\cr
}

The primary use of \pkg{secrdesign} is to predict by Monte Carlo
simulation the precision or bias of density estimates from different
detector layouts, given pilot values for density and the detection
parameters lambda0/g0 and sigma.

The important functions in \pkg{secrdesign} are:
\tabular{ll}{
  \code{\link{make.scenarios}} \tab generate dataframe of parameter values etc. \cr
  \code{\link{run.scenarios}} \tab perform simulations, with or without model fitting \cr
  \code{\link{fit.models}} \tab fit SECR model(s) to rawdata output from \code{run.scenarios} \cr
  \code{\link{predict.fittedmodels}} \tab infer `real' parameter
  estimates from fitted models \cr
  \code{\link{select.stats}} \tab collect output for a particular
  parameter \cr
  \code{\link{summary.selectedstatistics}} \tab numerical summary of results \cr 
  \code{\link{plot.selectedstatistics}} \tab histogram or CI plot for
  each scenario \cr 
}

Documentation is provided in a vignette
\url{../doc/secrdesign-vignette.pdf}; an Appendix has code for various
examples that should help get you started. The help pages are also
available as \url{../doc/secrdesign-manual.pdf}.

}

\seealso{

\code{\link{make.grid}}, \code{\link[secr]{sim.popn}}, \code{\link[secr]{sim.capthist}},  \code{\link[secr]{secr.fit}}
 
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\keyword{ package }
