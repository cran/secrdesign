\name{select.stats}
\alias{select.stats}
\alias{find.param}
\alias{find.stats}
  
\title{ Select Statistics to Summarize }

\description{

  When the results of each simulation with \code{run.scenarios} are
  saved as a dataframe (e.g. from predict()) it is necessary to select
  estimates of just one parameter for numerical summarization. This does
  the job. \code{find.param} is a helper function to quickly display the
  parameters available for summarisation.
  
}

\usage{

select.stats(object, parameter = "D", statistics, true)
find.param(object)
find.stats(object)

}

\arguments{
  \item{object}{`estimatetables' object from \code{\link{run.scenarios}}}
  \item{parameter}{character name of parameter to extract}
  \item{statistics}{character vector of statistic names}
  \item{true}{numeric vector of `'true' values of parameter, one per scenario}
}

\details{
  
\code{select.stats} is used to select a particular vector of numeric
values for summarization. The `parameter' argument indexes a row in the
data.frame for one replicate (i.e., one `real' parameter). Each
`statistic' is either a column in that data.frame or a statistic derived
from a column.

If \code{statistics} is not specified, the default is to use all numeric
columns in the input (i.e., c(`estimate', `SE.estimate', `lcl', `ucl')
for predict and c(`beta', `SE.beta', `lcl', `ucl') for coef).

\code{statistics} may include any of `estimate', `SE.estimate', `lcl',
`ucl', 'true', `RB', `RSE', `COV' and `ERR' (for outputtype `coef' use `beta' and
`SE.beta' instead of `estimate and `SE.estimate'). `true' refers to the
known parameter value used to generate the data.

The computed statistics are:
\tabular{llll}{
  Statistic \tab Name \tab Value \tab \cr
  RB \tab Relative bias \tab (estimate - true) / true \tab \cr
  RSE \tab Relative SE \tab SE.estimate / estimate \tab \cr
  ERR \tab Absolute deviation \tab abs(estimate - true) \tab \cr
  COV \tab Coverage \tab (estimate > lcl) & (estimate < ucl) \tab\cr
}

`RB', `COV' and `ERR' relate an estimate to the known (true) value of
the parameter in \code{object$scenarios}. They are computed only when a
model has been fitted without method = `none'.

`COV' remains binary (0/1) in the output from \code{select.stats}; the
result of interest is the mean of this statistic across replicates (see
\code{\link{summary.secrdesign}}). Similarly, `ERR' is used with field `rms'
in \code{\link{summary.secrdesign}} to compute the
root-mean-squared-error RMSE.

\code{find.param} and \code{find.stats} may be used to `peek' at objects
of class `estimatetables' and `selectedstatistics' respectively to
recall the available parameter estimates or `statistics'.

An attempt is made to extract \code{true} automatically if it is not provided. This does not always work (e.g. with extractfn \code{region.N}, region differing from the mask, and a heterogeneous density model). Check this by including ``true'' as a statistic to summarise (see Examples).

}

\value{
  
  For \code{select.stats}, an object with class
  c(`selectedstatistics',`secrdesign', `list') suitable for numerical
  summarization with \code{\link{summary.selectedstatistics}}. The value
  of `parameter' is stored as an attribute.

  For \code{find.param}, a character vector of the names of parameters
  with estimates in \code{object}.
  
}

\seealso{

  \code{\link{run.scenarios}}, \code{\link{validate}}
  
}

\examples{
## using nrepl = 2 just for checking
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid()
tmp1 <- run.scenarios(nrepl = 2, trapset = traps1, scenarios = scen1,
    fit = TRUE, extractfn = secr::trim)
tmp2 <- predict(tmp1)
tmp3 <- select.stats(tmp2, 'D', c('estimate','true','RB','RSE','COV'))
summary(tmp3)
}

\keyword{manip}