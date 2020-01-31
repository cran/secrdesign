\name{plot.optimalSpacing}
\alias{plot.optimalSpacing}
\alias{print.optimalSpacing}

\title{Plot and print methods for optimalSpacing object}

\description{

Plotsor print results from optimalSpacing.

}

\usage{

\method{plot}{optimalSpacing}(x, add = FALSE, plottype = c("RSE", "nrm"), ...)
\method{print}{optimalSpacing}(x, ...)

}

\arguments{
  \item{x}{object from \code{\link{optimalSpacing}}}
  \item{add}{logical; if TRUE will add to existing plot}
  \item{plottype}{character code}
  \item{\dots}{other arguments for plot, lines or points}
}

\details{

If \code{type = "RSE"} then RSE(D-hat) is plotted against R (relative detector spacing), otherwise the expected numbers of individuals, recaptures and movements are plotted against R.

The \dots argument may be used to pass other plotting arguments to override defaults:

\tabular{lll}{
Function \tab Arguments \tab Note \cr
plot \tab `xlab', `ylab', `xlim', `ylim', `las', `xaxs', `yaxs' \tab add = FALSE \cr
points \tab `col', `cex', `pch' \tab optimum and simulated RSE \cr
lines \tab `col', `lwd', `lty' \tab rule-of-thumb RSE \cr
}

The print method removes attributes before printing.

}

\value{None}

\seealso{\code{\link{optimalSpacing}}}

\keyword{ hplot }
