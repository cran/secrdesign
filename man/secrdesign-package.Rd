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
Version: \tab 2.1.2\cr
Date: \tab 2014-05-02\cr
License: \tab GNU General Public License Version 2 or later\cr
}

The primary use of \pkg{secrdesign} is to predict by Monte Carlo
simulation the precision or bias of density estimates from different
detector layouts, given pilot values for density and the detection
parameters g0 and sigma.

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

Documentation is provided in a vignette \url{../doc/secrdesign.pdf}; an
Appendix has code for various examples that should help get you started.
% \figure{workflow2.png}{workflow}
\if{html}{\figure{workflow2.png}{secrdesign workflow}}
\if{latex}{\figure{workflow2.png}{options: width=5.5in, height=3.5in}}

\code{\link{run.scenarios}} produces an object of class c(x,
`secrdesign', `list') where x depends on the arguments `fit' and
`extractfn'. The options for x (the primary class) are

\tabular{ll}{
  Class \tab Output from \code{extractfn()} for each replicate -- \cr
  \code{fittedmodels} \tab complete (or trimmed) fitted model from \code{\link[secr]{secr.fit}} \cr
  \code{estimatetables} \tab table of estimated parameters \cr
  \code{selectedstatistics} \tab vector of numerical values \cr
}

When output from each replicate is saved from
\code{\link{run.scenarios}} as a full or trimmed fitted model, further
processing with \code{\link{predict.fittedmodels}} or related functions
(\code{\link{coef.fittedmodels}}, \code{\link{derived.SL}}) will be
needed before summarization. The default (when \code{fit = TRUE}) is to
bypass this step by saving output from \code{run.scenarios} in
`estimatetables' form.

The help pages are also available as \url{../doc/secrdesign-manual.pdf}.

}

\seealso{
\code{\link{make.grid}},
\code{\link{make.scenarios}},
\code{\link{run.scenarios}},
\code{\link{select.stats}},
\code{\link{predict.fittedmodels}},
\code{\link{summary.selectedstatistics}},
\code{\link{plot.selectedstatistics}},
\code{\link[secr]{sim.popn}},
\code{\link[secr]{sim.capthist}}
\code{\link[secr]{secr.fit}}
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

%\references{
% Efford, M. G. (in prep.) Design of capture--recapture studies.
% }

\keyword{ package }
