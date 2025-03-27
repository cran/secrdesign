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
Version: \tab 2.9.3\cr
Date: \tab 2025-03-27\cr
License: \tab GNU General Public License Version 2 or later\cr
}

The primary use of \pkg{secrdesign} is to predict by Monte Carlo
simulation the precision or bias of density estimates from different
detector layouts, given pilot values for density and the detection
parameters lambda0/g0 and sigma.

Tools are also provided for predicting the performance of detector layouts 
without simulation, and for optimising layouts to meet various criteria, 
particularly expected counts.

The simulation functions in \pkg{secrdesign} are:
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

\if{html}{\figure{secrdesign-fig1.png}{options: width=700}}
\if{latex}{\figure{secrdesign-fig1.png}{options: width=5.5in}}

\bold{Fig.} Core simulation functions in \pkg{secrdesign} (yellow) and their main inputs and outputs. Output from the simulation function run.scenarios() may be saved as whole fitted models, predicted values (parameter estimates), or selected statistics. Each form of output requires different subsequent handling. The default path is shown by solid blue arrows. \cr

Other functions not used exclusively for simulation are:
\tabular{ll}{
  \code{\link{Enrm}} \tab expected numbers of individuals \eqn{n}, re-detections \eqn{r} and movements \eqn{m} \cr
  \code{\link{En2}} \tab expected number of individuals detected at two or more detectors \cr
  \code{\link{minnrRSE}} \tab approximate RSE(D-hat) given sample size (\eqn{n}, \eqn{r}) (Efford and Boulanger 2019)\cr
  \code{\link{GAoptim}} \tab optimization of detector placement using genetic algorithm (Durbach et al. 2021) \cr
  \code{\link{costing}} \tab various cost components \cr
  \code{\link{saturation}} \tab expected detector saturation (trap success) \cr
  \code{\link{scenarioSummary}} \tab applies \code{Enrm}, \code{minnrRSE}, and other summaries to each scenario in a dataframe \cr
  \code{\link{optimalSpacing}} \tab optimal detector spacing by rule-of-thumb and simulation RSE(D-hat)\cr
  %\code{\link{scenariosFromStatistics}} \tab match specified \eqn{n}, \eqn{r} \cr
}

A vignette documenting the simulation functions is available at \href{https://www.otago.ac.nz/density/pdfs/secrdesign-vignette.pdf}{secrdesign-vignette.pdf}. An Appendix in that vignette has code for various examples that should help get you started. 

Documentation for expected counts is in \href{https://www.otago.ac.nz/density/pdfs/secrdesign-Enrm.pdf}{secrdesign-Enrm.pdf}. Another vignette \href{https://www.otago.ac.nz/density/pdfs/secrdesign-tools.pdf}{secrdesign-tools.pdf} demonstrates other tools. These include the \code{optimalSpacing} function, for finding the detector spacing that yields the greatest precision for a given detector geometry, number of sampling occasions, density and detection parameters.

Help pages are also available as \url{../doc/secrdesign-manual.pdf}.

}

\references{

Durbach, I., Borchers, D., Sutherland, C. and Sharma, K. (2021) Fast, flexible 
alternatives to regular grid designs for spatial capture--recapture. 
\emph{Methods in Ecology and Evolution} \bold{12}, 298--310. \doi{10.1111/2041-210X.13517}

Efford, M. G., and Boulanger, J. (2019) Fast evaluation of study designs for 
spatially explicit capture--recapture. 
\emph{Methods in Ecology and Evolution}, \bold{10}, 1529--1535. \doi{10.1111/2041-210X.13239}

}

\seealso{

\code{\link[secr]{make.grid}}, 
\code{\link[secr]{sim.popn}}, 
\code{\link[secr]{sim.capthist}},  
\code{\link[secr]{secr.fit}}
 
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\keyword{ package }
