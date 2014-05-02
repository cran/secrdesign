\name{validate}
\alias{validate}

\title{Reject Implausible Statistics}

\description{
  
Simulation output may contain rogue values due to idiosyncracies of
model fitting. For example, nonidentifiability due to inadequate data
can result in spurious extreme `estimates' of the sampling
variance. Undue influence of rogue replicates can be reduced by using
the median as a summary field rather than the mean. This function is
another way to deal with the problem, by setting to NA selected
statistics from replicates for which some `test' statistic is
out-of-range.

}

\usage{
validate(x, test, validrange = c(0, Inf), targets = test)
}

\arguments{
  \item{x}{object that inherits from `selectedstatistics'}
  \item{test}{character; name of statistic to check}
  \item{validrange}{numeric vector comprising the minimum and maximum
  permitted values of `test', or a matrix (see details)}
  \item{targets}{character vector with names of one or more statistics
  to set to missing (NA) when test is out-of-range}
}

\details{

  Values of `test' and `targets' should be columns in each component
  `replicate x statistic' matrix (i.e., scenario) of
  \code{x$output}. You can check for these with
  \code{\link{find.stats}}.

  If \code{validrange} is a matrix its first and second columns are
  interpreted as scenario-specific bounds (minima and maxima), and the
  number of rows must match the number of scenarios.
  
If all non-missing values of `test' are in the valid range, the effect
is to force the target statistics to NA wherever `test' is NA.

The default is to change only the test field itself. If the value of
`test' does not appear in `targets' then the test field is unchanged.

If \code{targets = "all"} then all columns are set to NA when the test
fails.

}

\value{

  An object of class c(`selectedstatistics', secrdesign', `list') with
  the same structure and header information as the input, but possibly
  with some values in the `output' component converted to NA.

}

\seealso{
  
  \code{\link{select.stats}}, \code{\link{find.stats}}
  
}

\examples{

\dontrun{

## generate some data
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid()
tmp1 <- run.scenarios(nrepl = 5, trapset = traps1, scenarios = scen1,
    fit = TRUE, extractfn = trim)
tmp2 <- predict(tmp1)
tmp3 <- select.stats(tmp2, 'D', c('estimate','RB','RSE','COV'))

## just for demonstration --
## apply scenario-specific +/- 20\% bounds for estimated density
## set RB, RSE and COV to NA when estimate is outside this range
permitted <- outer(tmp3$scenarios$D, c(0.8,1.2)) 
permitted   ## a 2 x 2 matrix
tmp4 <- validate(tmp3, 'estimate', permitted, c('RB', 'RSE','COV'))

## what have we done?!
tmp4$output
summary(tmp4)

}
}

\keyword{ manip }

