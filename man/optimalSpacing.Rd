\name{optimalSpacing}
\alias{optimalSpacing}

\title{Optimal Detector Spacing}

\description{Estimate the detector spacing that yields the greatest precision for a given detector geometry, number of sampling occasions, density and detection parameters.}

\usage{

optimalSpacing (D, traps, detectpar, noccasions, nrepeats = 1, 
    detectfn = c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX'),
    fittedmodel = NULL, xsigma = 4, R = seq(0.2, 4, 0.2), CF = 1.0,
    distribution = c("poisson", "binomial"),
    fit.function = c("none", "secr.fit"), 
    simulationR = seq(0.4, 4, 0.4), nrepl = 10, 
    plt = FALSE, ...) 
}

\arguments{

  \item{D}{population density animals / hectare (constant)}
  \item{traps}{ \code{\link{traps}} object }
  \item{detectpar}{named list giving a value for each parameter of detection function (sigma not needed)}
  \item{noccasions}{integer number of sampling occasions}
  \item{nrepeats}{integer number of replicate arrays (not yet used)}
  \item{detectfn}{integer code or character string for shape of detection function -- see \link{detectfn}}
  \item{fittedmodel}{secr fitted model (instead of preceding arguments)}
  \item{xsigma}{numeric buffer width as multiple of sigma}
  \item{R}{numeric vector of relative spacings at which to plot rule-of-thumb RSE(D-hat)}
  \item{CF}{numeric correction factor for rule-of-thumb RSE}
  \item{distribution}{character distribution of number of individuals detected}
  \item{fit.function}{character function to use for model fitting}
  \item{simulationR}{numeric vector of relative spacings at which to simulate}
  \item{nrepl}{integer number of replicate simulations }
  \item{plt}{logical; if TRUE then results are plotted}
  \item{\dots}{other arguments passed to various functions (see Details)}
  
}

\details{

A numerical search over possible spacings uses the rule-of-thumb RSE(D-hat) given by \code{\link{minnrRSE}} as the objective function. 

\code{traps} provides the geometry of the detector layout and the initial spacing \eqn{s}. Function \code{\link{optimize}} is used to search for a solution (minimum RSE) in the range of \code{R} x \eqn{s}.

The computation emulates variation in detector spacing by inverse variation in sigma (sigma' = sigma / R) with compensating variation in density. Mask buffer width and spacing are also scaled by R.

If \code{fit.function} is not "none" then simulations are also performed for the relative spacings in \code{simulationR}. Density, sigma and mask attributes are scaled as for the rule-of-thumb calculations. Using `method = "none"` gives fast prediction of RSE (from the Hessian evaluated at the known parameter values), but does not estimate bias.

The \dots argument may be used to set the values of these arguments:

\tabular{ll}{
Function \tab Arguments \cr
make.mask \tab `nx', `type', `poly',`poly.habitat' \cr
run.scenarios \tab `seed', `ncores', `method' \cr
plot.optimalSpacing \tab `add', \dots \cr
}

The argument \code{CF} may be set to NA to suppress rule-of-thumb RSE, including optimisation. \code{range(R)} specifies the search interval for optimisation.

A \code{plot} method is provided, with options for plotting different components.

}

\value{

List of two components, one for the rule-of-thumb optimisation (rotRSE) and the other for simulation results, if requested (simRSE).

The optimisation results are
\item{values}{dataframe with E(n), E(r) and the rule-of-thumb RSE for each requested R}
\item{optimum.spacing}{the absolute spacing that yields maximum precision (minimum rule-of-thumb RSE(D-hat))}
\item{optimum.R}{spacing relative to sigma}
\item{minimum.RSE}{final value of the objective function (minimum rule-of-thumb RSE(D-hat))}

The simulation results in the dataframe simRSE are the mean and SE of the simulated RSE(D-hat) for each level of \code{simulationR}, with added columns for the relative bias (RB) and relative root-mean-square-error (rRMSE) of D-hat.

Results are returned invisibly if \code{plt = TRUE}.

}

\note{fit.function = 'openCR.fit' was deprecated from 2.5.8 and has been removed as an option}

\seealso{

\code{\link{minnrRSE}},
\code{\link{plot.optimalSpacing}}

}

\section{Warnings}{

For single-catch traps, use of a maximum likelihood estimate of lambda0 from a fitted multi-catch model results in negative bias.

Only hazard-based detection functions are supported. The meaning of the `sigma' parameter depends on the function, and so will the optimal spacing in sigma units.
}

\examples{

grid <- make.grid(7, 7)  # default multi-catch detector
optimalSpacing(D = 5, traps = grid, detectpar = list(lambda0 = 0.2, sigma = 20), 
   noccasions = 5, plt = TRUE)

\dontrun{

optimalSpacing(D = 5, traps = grid, detectpar = list(lambda0 = 0.4, sigma = 20), 
    detectfn = 'HEX', R = seq(1,6,0.4), noccasions = 10, plt = TRUE, col = "blue")

## with simulations
grid <- make.grid(8, 8, spacing = 20, detector = 'proximity')
optimalSpacing(D = 5, traps = grid, detectfn = "HHN", detectpar = 
    list(lambda0 = 0.2, sigma = 20), noccasions = 5, nrepl = 20, nx = 32, 
    ncores = 4, plt = TRUE, col = "blue")

## manual check
grid <- make.grid(8, 8, spacing = 60, detector = 'proximity')
scen <- make.scenarios(D = 5, detectfn = 14, lambda0 = 0.2, sigma = 20,
    noccasions = 5)
sim1 <- run.scenarios(nrepl = 20, scen, trapset = list(grid), fit = TRUE, 
    fit.args = list(detectfn = 14), ncores = 4, byscenario = FALSE)
summary(sim1)

}

}


\keyword{ design }

