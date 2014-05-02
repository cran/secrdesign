\name{make.array}
\alias{make.array}
  
\title{ Re-cast Simulated Statistical Output as Array }

\description{

This function is used internally by \code{\link{summary.secrdesign}}, and may
occasionally be of general use.  
}

\usage{

make.array(object)

}

\arguments{ \item{object}{secrdesign object containing numerical values
  for a particular parameter (i.e. output from
  \code{\link{select.stats}} 
  inheriting from `selectedstatistics') } }

\details{
  
\code{make.array} converts a particular simulated numerical output into an array
with one dimension for each varying input.

}

\value{
  
  A numeric array with dimensions corresponding to the varying inputs.
  
}

\seealso{

  \code{\link{run.scenarios}}
  
}

\examples{
## collect raw counts
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid()
tmp1 <- run.scenarios(nrepl = 50, trapset = traps1, scenarios = scen1,
    fit = FALSE)
make.array(tmp1)
}

\keyword{manip}