\name{predict.fittedmodels}
\alias{predict.fittedmodels}
\alias{coef.fittedmodels}
\alias{derived.SL}
\alias{regionN.SL}
  
\title{ Extract Estimates From Fitted Models }

\description{

  If simulations have been saved from \code{run.scenarios} as fitted
  secr models it is necessary to use one of these functions to extract
  estimates for later summarization.
  
}

\usage{

\method{predict}{fittedmodels}(object, ...)

\method{coef}{fittedmodels}(object, ...)

derived.SL(object, ...)

regionN.SL(object, ...)

}

\arguments{
  \item{object}{fitted model simulation output from \code{\link{run.scenarios}}}
  \item{\dots}{other arguments passed to predict, coef, derived or region.N }
}

\details{
  
These functions are used when output from \code{\link{run.scenarios}}
has been saved as fitted models. \code{derived.SL} and \code{regionN.SL}
require a full fit (including the design object) whereas a trimmed model
is sufficient for \code{predict} and \code{coef}.

\code{derived.SL} is used to compute the Horvitz-Thompson-like estimate
of density when \code{\link[secr]{secr.fit}} has been used with \code{CL
= TRUE}; it is roughly equivalent to \code{predict}.

\code{regionN.SL} predicts the realised number (R.N) or expected number
(E.N) in a masked area. When detector layouts and/or \code{sigma} vary,
the masked area will also vary (arbitrarily, depending on the buffer
argument `xsigma') unless a mask is provided by the user; this may be
done either in \code{run.scenarios} or in \code{regionN.SL}.

}

\value{
  
  An object with class (`estimatetables', `secrdesign', `list') with
  appropriate outputtype (`predicted', `coef', `derived', `regionN';
  sSee also \code{\link{run.scenarios}}).

}

\seealso{

  \code{\link{run.scenarios}}
  
}

\examples{
## using nrepl = 2 just for checking
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid()  ## default 6 x 6 grid of multi-catch traps
tmp1 <- run.scenarios(nrepl = 2, trapset = traps1, scenarios = scen1,
    fit = TRUE, extractfn = trim)
tmp2 <- predict(tmp1)
tmp3 <- select.stats(tmp2, 'D', c('estimate','RB','RSE'))
summary(tmp3)
}

\keyword{manip}