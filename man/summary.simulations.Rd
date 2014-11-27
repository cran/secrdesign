\name{summary.secrdesign}
\alias{summary.secrdesign}
\alias{summary.rawdata}
\alias{summary.estimatetables}
\alias{summary.selectedstatistics}
\alias{plot.selectedstatistics}
\alias{header}

\title{Generic Methods for secrdesign Objects}

\description{
  
  Methods to summarize simulated datasets.
  
}

\usage{

\method{summary}{secrdesign}(object, ...)

\method{summary}{rawdata}(object, ...)

\method{summary}{estimatetables}(object, ...)

\method{summary}{selectedstatistics}(object, fields = c('n', 'mean',
'se'), dec = 5, alpha = 0.05, type = c('list','dataframe','array'), ...)

\method{plot}{selectedstatistics}(x, scenarios, statistic, type =
c('hist', 'CI'), refline, xlab = NULL, ...)

header(object)

}

\arguments{
  \item{object}{object of class simulations from \code{run.scenarios}}
  \item{dec}{number of decimal places in output}
  \item{fields}{character vector; names of required summary
    statistics (see Details)}
  \item{alpha}{ alpha level for confidence intervals and quantiles}
  \item{type}{character code for type of output (see Details)}
  \item{\dots}{other arguments -- not currently used by summary but
  passed to \code{\link{hist}} by the plot method}
  \item{x}{object of class `selectedstatistics' from
    \code{run.scenarios}}
  \item{scenarios}{integer indices of scenarios to plot (all plotted if not specified)}
  \item{statistic}{integer or character indices if the statistics in x for which
    histograms are requested}
  \item{refline}{logical; if TRUE a reference line is plotted at the
  true value of a parameter}
  \item{xlab}{character; optional label for x-axis}

}

\details{

    If \code{object} inherits from `selectedstatistics' then the numeric

  If \code{object} inherits from `selectedstatistics' then the numeric
  results from replicate simulations are summarized using the chosen
  `fields' (by default, the number of non-missing values, mean and standard
  error), along with header information describing the
  simulations. Otherwise the header alone is returned.
  
  \code{fields} is a vector of any selection from c(`n', `mean', `sd',
  `se', `min', `max', `lcl', `ucl', `median', `q', `rms'), or the
  character value `all'.

  Field `q' provides 1000 \code{alpha/2} and 1000[1 - \code{alpha/2}]
  quantiles qxxx and qyyy.

  `lcl' and `ucl' refer to the upper and lower limits of a 100(1 -
  alpha)\% confidence interval for the statistic, across replicates.

  `rms' gives the root-mean-square of the statistic - most useful for
  the statistic `ERR' (see \code{\link{select.stats}}) when it
  represents the overall accuracy or RMSE.

  The \code{plot} method plots either (i) histograms of the selected
  statistics (type = `hist') or (ii) the estimate and confidence interval for
  each replicate (type = `CI'). The default for type = `hist' is to plot
  the first statistic - this is usually `n' (number of detected animals)
  when \code{fit = FALSE}, and `estimate' (parameter estimate) when
  \code{fit = TRUE}. If length(statistic) > 1 then more than one plot
  will be produced, so a multi-column or multi-row layout should be
  prepared with \code{par} arguments `mfcol' or `mfrow'.

  For type = `CI' the statistics must include `estimate', `lcl' and
  `ucl' (or `beta', `lcl' and `ucl' if outputtype = `coef').
  
}

\value{
  List with components `header'
 
  \item{call}{original function call}
  \item{starttime}{from object}
  \item{proctime}{from object}
  \item{constants}{small dataframe with values of non-varying inputs}
  \item{varying}{small dataframe with values of varying inputs}
  \item{fit.args}{small dataframe with values arguments for
  secr.fit, if specified}

  and `OUTPUT', a list with one component for each field. Each component
  may be a list or an array.
  
}

\seealso{
  
  \code{\link{run.scenarios}},
  \code{\link{make.array}},
  \code{\link{select.stats}}
  \code{\link{validate}}

}

\examples{

## collect raw counts
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid()
tmp1 <- run.scenarios(nrepl = 50, trapset = traps1, scenarios = scen1,
    fit = FALSE)

opar <- par(mfrow=c(2,3))
plot(tmp1, statistic = 1:3)
par(opar)

summary(tmp1)

summary(tmp1, field=c('q025', 'median', 'q975'))

}
\keyword{ Generic }