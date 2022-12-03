\name{Internal}

\alias{compactSample}
\alias{GApenfn}

\title{Internal Functions}

\description{ 

Functions that are called internally by \pkg{secrdesign}. 
These are exported and may be called separately for testing.

}

\usage{

compactSample (traps, n)

GApenfn(traps, sigma)

}

\arguments{

  \item{traps}{secr trapsobject}
  \item{n}{integer number in sample (\eqn{0 < n \le ntraps)}}
  \item{sigma}{numeric sparial scale parameter}
  
}

\details{

  \code{compactSample} selects a detector at random and returns the a compact subset of surrounding detectors.

  \code{GApenfn} is the default pen_fn used by \code{\link{GAoptim}} When called with a non-null penalty argument.
  
}

\value{

GApenfn -- a numeric vector with the number of trap pairs separated by 2.5-3.5 sigma and 3.5-4.5 sigma.

compactSample -- an object like traps, but with only n rows.
}

\seealso{

\code{\link{GAoptim}}, 

}

\references{

Durbach, I., Borchers, D., Sutherland, C. and Sharma, K. (2021) Fast, flexible alternatives to regular grid designs for spatial capture--recapture. \emph{Methods in Ecology and Evolution} \bold{12}, 298--310. DOI 10.1111/2041-210X.13517

}

\examples{

CStraps <- compactSample(traps(captdata), n = 20)

plot(traps(captdata))
plot(CStraps, add = TRUE, detpar = list(fg = 'blue',pch = 16))

GApenfn(CStraps, sigma = 25)

}