\name{minsimRSE}
\alias{minsimRSE}
\alias{minsimRSE.optimalSpacing}

\title{ Optimal Spacing by Simulation }

\description{

A method to obtain a unique `optimal' spacing from previously simulated scenaios for detector spacing.

}

\usage{

\method{minsimRSE}{optimalSpacing}(object, cut = 0.2, plt = FALSE, verbose = FALSE, incr = 0.1, ...)

}

\arguments{

  \item{object}{optimalSpacing object}
  \item{cut}{numeric maximum \eqn{\Delta} RSE to include}
  \item{plt}{logical; if TRUE a plot is generated}
  \item{verbose}{logical; if TRUE then output includes fitted model}
  \item{incr}{numeric spacing of computed points (R)}
  \item{\dots}{other arguments passed to \code{plot.optimalSpacing}}

}

\details{

A quadratic is fitted to the simulated RSE (y) vs simulationR (x), including only values of x and y for which \eqn{y \le \mbox{min}(y) \times (1+\mbox{cut})}. The restriction allows the user to exclude extreme x-values for which the quadratic is a poor fit.

The optimum is the minimum of the quadratic \eqn{ax^2 + bx + c}, given by \eqn{-b/2a}.

The quadratic is fitted with \code{\link{lm}} (\code{lm(RSE.mean ~ R + I(R^2)}).

}

\value{

When \code{verbose = FALSE}, a numeric vector with optimum R (multiple of sigma) and corresponding RSE.

When \code{verbose = TRUE}, a list with components --

  \item{model}{fitted model from \code{\link{lm}}}
  \item{fitted}{dataframe of points on fitted curve}
  \item{R}{optimum R}
  \item{RSE}{minimum RSE}

}

\seealso{

\code{\link{optimalSpacing}}

}

\examples{

\donttest{

grid <- make.grid(8, 8, spacing = 20, detector = 'proximity')

# method = "none" uses the shortcut variance
tmp <- optimalSpacing(D = 5, traps = grid, detectfn = "HHN", detectpar = 
    list(lambda0 = 1, sigma = 20), noccasions = 1, nx = 32, 
    fit.function = "secr.fit", method = "none", simulationR = seq(1.2,2.2,0.2))
minsimRSE(tmp, plt = TRUE)

}

}
