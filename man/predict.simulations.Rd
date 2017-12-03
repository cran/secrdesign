\name{predict.fittedmodels}
\alias{predict}
\alias{coef}
\alias{derived}
\alias{region.N}
\alias{predict.fittedmodels}
\alias{coef.fittedmodels}
\alias{derived.fittedmodels}
\alias{region.N.fittedmodels}
  
\title{ Extract Estimates From Fitted Models }

\description{

  If simulations have been saved from \code{run.scenarios} as fitted
  secr models it is necessary to use one of these functions to extract
  estimates for later summarization.
  
}

\usage{

\method{predict}{fittedmodels}(object, ...)

\method{coef}{fittedmodels}(object, ...)

\method{derived}{fittedmodels}(object, ...)

\method{region.N}{fittedmodels}(object, ...)

}

\arguments{
  \item{object}{fitted model simulation output from \code{\link{run.scenarios}}}
  \item{\dots}{other arguments passed to predict, coef, derived or region.N }
}

\details{
  
These functions are used when output from \code{\link{run.scenarios}}
has been saved as fitted models. \code{derived} and \code{region.N}
require a full fit (including the mask and design0 objects) whereas a trimmed model
is sufficient for \code{predict} and \code{coef}.

\code{derived} is used to compute the Horvitz-Thompson-like estimate
of density when \code{\link[secr]{secr.fit}} has been used with \code{CL
= TRUE}; it is roughly equivalent to \code{predict}.

\code{region.N} predicts the realised number (R.N) or expected number
(E.N) in a masked area. When detector layouts and/or \code{sigma} vary,
the masked area will also vary (arbitrarily, depending on the buffer
argument `xsigma') unless a mask is provided by the user; this may be
done either in \code{run.scenarios} or in \code{region.N}.

}

\value{
  
  An object with class (`estimatetables', `secrdesign', `list') with
  appropriate outputtype (`predicted', `coef', `derived', `regionN';
  see also \code{\link{run.scenarios}}).

}

\seealso{

  \code{\link{run.scenarios}}
  \code{\link[secr]{coef.secr}}
  \code{\link[secr]{predict.secr}}
  \code{\link[secr]{derived.secr}}
  \code{\link[secr]{region.N.secr}}
}

\examples{

\dontrun{
scen1 <- make.scenarios(D = c(3,6), sigma = 25, g0 = 0.2)
traps1 <- make.grid()  ## default 6 x 6 grid of multi-catch traps
tmp1 <- run.scenarios(nrepl = 10, trapset = traps1, scenarios = scen1,
    fit = TRUE, extractfn = trim)
tmp2 <- predict(tmp1)
tmp3 <- select.stats(tmp2, 'D', c('estimate','RB','RSE'))
summary(tmp3)

## for derived and region.N need more than just 'trimmed' secr object
## use argument 'keep' to save mask and design0 usually discarded by trim
tmp4 <- run.scenarios(nrepl = 10, trapset = traps1, scenarios = scen1,
    fit = TRUE, extractfn = trim, keep = c('mask','design0'))

summary(derived(tmp4))

## for region.N we must specify the parameter for which we want statistics
## (default 'D' not relevant)
tmp5 <- select.stats(region.N(tmp4), parameter = 'E.N')
summary(tmp5)
}

}

\section{Note}{

From \pkg{secrdesign} 2.5.3 the methods described here replace the functions \code{derived.SL} and \code{regionN.SL}. This is for compatibility with \pkg{secr}.

}

\keyword{manip}