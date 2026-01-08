\name{plot.optimalSpacing}
\alias{plot.optimalSpacing}
\alias{print.optimalSpacing}

\title{Plot and print methods for optimalSpacing object}

\description{

Plotsor print results from optimalSpacing.

}

\usage{

\method{plot}{optimalSpacing}(x, add = FALSE, plottype = c("both", "RSE", "nrm", "RB", "RMSE"), 
    xtype = c('relative','absolute'), xoffset = 0, ...)
\method{print}{optimalSpacing}(x, ...)

}

\arguments{
  \item{x}{object from \code{\link{optimalSpacing}}}
  \item{add}{logical; if TRUE will add to existing plot}
  \item{plottype}{character code}
  \item{xtype}{character scale of x-axis}
  \item{xoffset}{numeric offset of points on x-axis}
  \item{\dots}{other arguments for plot, lines or points}
}

\details{

If \code{plottype = "RSE"} then RSE(D-hat) is plotted against R (relative detector spacing). \cr
If \code{plottype = "nrm"} then the expected numbers of individuals, recaptures and movements are plotted against R. \cr
If \code{plottype = "RB"} then the simulated relative bias of D-hat is plotted against R. \cr
If \code{plottype = "rRMSE"} then the simulated relative root-mean-square-error of D-hat is plotted against R.

\code{plottype = "both"} then both the simulated and approximate RSE are plotted.

If \code{xtype = "absolute"} then the scale of the x-axis is R * sigma (from the detectpar attribute of x).

\code{xoffset} is used to slightly displace added points to avoid overlap.

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
