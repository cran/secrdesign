\name{transformOutput}
\alias{transformOutput}

\title{

Transform Simulation Output

}

\description{

Transform output component of simulation output from \code{\link{run.scenarios}}. 
Typically this replaces an entire saved model fit with a table of estimates from that fit.

}

\usage{

transformOutput(object, extractfn, outputtype = "predicted", ...)

}

\arguments{
  \item{object}{output from \code{\link{run.scenarios}}}
  \item{extractfn}{function such as the `extractfn' argument of \code{\link{run.scenarios}}}
  \item{outputtype}{character (see Details)}
  \item{\dots}{other arguments passed to \code{extractfn}}
}

\details{

Each replicate of each scenario is transformed using `extractfn', which should accept as input the object returned by the extractfn of the original call to \code{run.scenarios}. As a typical example, \code{sims <- run.scenarios(..., fit = TRUE, extractfn = identity)} returns outputs of class `secr' and could be followed by \code{sims2 <- transformOutput(sims, predict)}; \code{sims2} may be used as input to \code{\link{estimateSummary}} and other summary functions.

}

\value{

An object resembling the output from \code{\link{run.scenarios}} but with transformed output.
The outputtype and class of the object are changed to match `outputtype'.

}

\seealso{

\code{\link{run.scenarios}},
\code{\link{estimateSummary}},
\code{\link{outputtype<-}}

}
