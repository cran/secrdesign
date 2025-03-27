\name{Internal}

\alias{compactSample}
\alias{GApenfn}
\alias{outputtype}
\alias{outputtype<-}

\title{Internal Functions}

\description{ 

Functions that are called internally by \pkg{secrdesign}. 
These are exported and may be called separately for testing etc.

}

\usage{

compactSample (traps, n)

GApenfn(traps, sigma)

'outputtype<-'(object, value)

}

\arguments{

  \item{traps}{secr trapsobject}
  \item{n}{integer number in sample (\eqn{0 < n \le ntraps)}}
  \item{sigma}{numeric sparial scale parameter}
  \item{object}{object output from \code{\link{run.scenarios}}}
  \item{value}{replacement value for outputtype of object}
}

\details{

  \code{compactSample} selects a detector at random and returns the a compact subset of surrounding detectors.

  \code{GApenfn} is the default pen_fn used by \code{\link{GAoptim}} When called with a non-null penalty argument.
  
  Values of outputtype map to class of the run.scenarios output as follows
  \tabular{lll}{
  Output type \tab Class \cr
  secrfit \tab c("fittedmodels", "secrdesign", "list") \cr
  ipsecrfit \tab c("fittedmodels", "secrdesign", "list") \cr
  predicted \tab c("estimatetables", "secrdesign", "list") \cr
  derived \tab c("estimatetables", "secrdesign", "list") \cr
  regionN \tab c("estimatetables", "secrdesign", "list") \cr
  coef \tab c("estimatetables", "secrdesign", "list") \cr
  user \tab c("estimatetables", "secrdesign", "list") \cr
  secrsummary \tab c("summary", "secrdesign", "list") \cr
  capthist \tab c("rawdata", "secrdesign", "list") \cr
  selectedstatistics \tab c("selectedstatistics", "secrdesign", "list") \cr
  }
  
Calling the replacement function automatically changes the class of the output object as appropriate. This determines how the output is handled by downstream functions such as \code{summary}. Using a custom extractfn or post-processing the output sometimes requires the outputtype to be set manually (see example in the Multi-model section of secrdesign-vignette.pdf).
  
}

\value{

GApenfn -- a numeric vector with the number of trap pairs separated by 2.5-3.5 sigma and 3.5-4.5 sigma.

compactSample -- an object like traps, but with only n rows.
}

\seealso{

\code{\link{GAoptim}}, 

}

\references{

Durbach, I., Borchers, D., Sutherland, C. and Sharma, K. (2021) Fast, flexible 
alternatives to regular grid designs for spatial capture--recapture. 
\emph{Methods in Ecology and Evolution} \bold{12}, 298--310. \doi{10.1111/2041-210X.13517}

}

\examples{

CStraps <- compactSample(traps(captdata), n = 20)

plot(traps(captdata))
plot(CStraps, add = TRUE, detpar = list(fg = 'blue',pch = 16))

GApenfn(CStraps, sigma = 25)

}