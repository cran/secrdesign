\name{rbind.estimatetables}
\alias{rbind.estimatetables}
\alias{rbind.selectedstatistics}
\alias{c.estimatetables}
\alias{c.selectedstatistics}

\title{

Combine Simulation Output

}

\description{

Methods to combine output from separate executions of \code{\link{run.scenarios}}.

}

\usage{

\method{rbind}{estimatetables}(..., deparse.level = 1)

\method{rbind}{selectedstatistics}(..., deparse.level = 1)

\method{c}{estimatetables}(...)

\method{c}{selectedstatistics}(...)

}

\arguments{
  \item{\dots}{estimatetables or selectedstatistics output from \code{\link{run.scenarios}}}
  \item{deparse.level}{not used (required by generic method \code{\link{rbind}})}
}

\details{
\code{rbind} assumes all inputs used exactly the same scenarios. Replicate estimate tables are combined across executions for each scenario in turn. This is useful to increase the number of replicates by combining two batches of simulations with different random seeds. The `scenarios' component remains unchanged.

\code{c} combines outputs from \code{run.scenarios} that may differ in their scenarios. The `output' component of the result is a concatenation of the output lists in the input. The `scenarios' component of the result comprises the input scenarios stacked with \code{\link{rbind.data.frame}}.

The compatibility of the inputs is checked, but the checks are not exhaustive. Users should be wary.

}

\value{

`estimatetables' or `selectedstatistics' object combining the inputs

}

\seealso{

\code{\link{make.scenarios}}
\code{\link{run.scenarios}}

}

\examples{

## Simple example: generate and summarise trapping data at two densities
## result inherits from 'selectedstatistics'

scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2, noccasions = 5)
traps1 <- make.grid()   ## default 6 x 6 trap grid
tmp1 <- run.scenarios(nrepl = 5, trapset = traps1, scenarios = scen1,
    fit = FALSE, seed = 123)
tmp2 <- run.scenarios(nrepl = 15, trapset = traps1, scenarios = scen1,
    fit = FALSE, seed = 127)
    
summary(rbind(tmp1,tmp2))    
summary(c(tmp1,tmp2))    

}
