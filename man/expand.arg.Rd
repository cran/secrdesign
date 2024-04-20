\name{expand.arg}
\alias{expand.arg}

\title{Expand Argument List}

\description{

Generate a list of lists from vectors of argument values.

}

\usage{

expand.arg(..., sublist = list())

}

\arguments{

  \item{\dots}{named vectors of argument values}
  \item{sublist}{named list of character vectors}

}

\details{

The full 'pop.args', 'det.args' and 'fit.args' arguments of \code{\link{run.scenarios}} are lists of lists corresponding to the popindex, detindex and fitindex columns in scenarios.

\code{expand.arg} constructs such lists from all possible combinations of specified arguments, with invariant arguments appended. 

\code{sublist} may be specified to shift one or more named subsets of \dots arguments to a sub-list such as 'detectpar' (see Examples).

}

\value{

A list of lists. A flat dataframe of combinations is returned as the attribute 'comb'.

}

\seealso{

\code{\link{run.scenarios}}

}

\examples{

## generate some arguments for sim.capthist
## noccasions is constant

expand.arg(detectfn = c('HN','EX'), noccasions = 5)

## detectpar sub-list

expand.arg(lambda0 = c(0.1, 0.2), sigma = 20, detectfn = c('HHN','HEX'),
    sublist = list(detectpar = c('lambda0','sigma')))

}

