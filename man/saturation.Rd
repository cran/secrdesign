\name{saturation}
\alias{saturation}

\title{Detector saturation}

\description{

Computes the expected proportion of successful detectors (i.e., `trap success'). The calculation does not allow for local variation in realised density (number of animals centred near each detector) and the predictions are therefore slightly higher than simulations with Poisson local density. The discrepancy is typically less than 1\%.

}

\usage{

saturation(traps, mask, detectpar, detectfn = 
    c("HHN", "HHR", "HEX", "HAN", "HCG", 'HN', 'HR', 'EX'), 
    D, plt = FALSE, add = FALSE, ...)

}

\arguments{

  \item{traps}{secr traps object}
  \item{mask}{secr mask object}
  \item{detectpar}{ a named list giving a value for each parameter of detection function }
  \item{detectfn}{integer code or character string for shape of detection function 
      -- see \link[secr]{detectfn}}
  \item{D}{population density animals / hectare; may be scalar or vector of length \code{nrow(mask)}}
  \item{plt}{logical; if TRUE then a colour plot is produced}
  \item{add}{logical; if TRUE any plot is added to the existing plot}
  \item{\dots}{other arguments passed to plot.mask when plt = TRUE}
  
}
\details{

The calculation is based on an additive hazard model. If \code{detectfn} is not a hazard function (`HHN', `HEX', `HHR', `HAN' and `HCG') then an attempt is made to approximate one of the hazard functions (HN -> HHN, HR -> HHR, EX -> HEX). The default is `HHN'.

Computation is not possible for single-catch traps.

An empirical estimate of saturation is the total number of detectors visited divided by the total number of detectors used. These are outputs from the summary method for capthist objects. See Examples.

}

\value{

A list with components
\item{bydetector}{expected saturation for each detector}
\item{mean}{average over detectors}

The list is returned invisibly if \code{plt = TRUE}.

}

\seealso{

\code{\link{Enrm}}

}

\examples{

tr <- traps(captdata)
detector(tr) <- 'multi'
mask <- make.mask(tr, buffer = 100)
saturation(tr, mask, detectpar = list(lambda0 = 0.27, sigma = 29), 
    detectfn = 'HHN', D = 5.5, plt = TRUE)
plotMaskEdge(as.mask(tr), add = TRUE)  ## boundary line

# empirical - useful for extractfn argument of secrdesign::run.scenarios
satfn <- function(CH) { 
    sumCH <- summary(CH)$counts
    sumCH['detectors visited', 'Total'] /  sumCH['detectors used', 'Total']
}
satfn(captdata)

}

\keyword{ manip }