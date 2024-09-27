\name{count}
\alias{count}
\alias{coef.summary}
\alias{predict.summary}
\alias{count.summary}

\title{Extract Summaries}

\description{

Reshape results from \code{run.scenarios(..., extractfn = summary)} so that they may be passed to the usual summary functions of \pkg{secrdesign}.

}

\usage{

count(object, ...)

\method{predict}{summary}(object, ...)
\method{coef}{summary}(object, ...)
\method{count}{summary}(object, ...)

}

\arguments{

  \item{object}{summary simulation output from \code{\link{run.scenarios}}}
  \item{\dots}{other arguments (not used)}

}

\details{

The aim is to extract numerical results from simulations performed using \code{run.scenarios(..., extractfn = summary)}. The results may then be passed to the summary method for `secrdesign' objects, possibly via \code{\link{select.stats}} (see Examples).

}

\value{

An object of class \code{c("estimatetables", "secrdesign", "list")} in which the output component for each scenario is a list of dataframes, one per replicate. The structure of each dataframe is indicated in the following table (parameters may vary with model); `parameters' and `statistics' correspond to arguments of \code{\link{select.stats}}.

\tabular{lllll}{
Function \tab \tab Row(s) \tab \tab Columns \cr
 \tab \tab (parameters) \tab \tab (statistics) \cr

count \tab \tab  Number \tab \tab Animals, Detections, Moves \cr
coef \tab \tab  D, g0, sigma \tab \tab  estimate, SE.estimate, lcl, ucl \cr
predict \tab \tab  D, g0, sigma \tab \tab  estimate, SE.estimate, lcl, ucl \cr

}

}

\seealso{

\code{\link[secr]{predict.secr}},
\code{\link[secr]{coef.secr}},

}

\examples{

## generate some simulations
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid(6, 6, spacing = 25)
sims1 <- run.scenarios(nrepl = 2, trapset = traps1, scenarios =
     scen1, seed = 345, fit = TRUE, extractfn = summary)

## view the results
count(sims1)$output
predict(sims1)$output

summary(sims1)  ## header only

summary(count(sims1))   # equivalent to following
summary(select.stats(count(sims1), parameter = 'Number'))  

summary(predict(sims1)) # default select.stats parameter = 'D'
summary(select.stats(predict(sims1), parameter = 'sigma')  )

}

\keyword{ datagen }
