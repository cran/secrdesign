\name{getdetectpar}
\alias{getdetectpar}

\title{Ballpark Detection Parameters}

\description{Detection parameters for an animal population may be guessed from some basic inputs (population density, a coefficent of home-range overlap, and the expected number of detections on a given detector array). These values are useful as a starting point for study design. They are not 'estimates'.
}

\usage{

getdetectpar(D, C, sigma = NULL, k = 0.5, ...)

}

\arguments{

  \item{D}{population density animals / hectare; may be scalar or vector of length \code{nrow(mask)}}
  \item{C}{integer expected total number of detections}
  \item{sigma}{numeric spatial scale parameter of chosen detection function, in metres (optional)}
  \item{k}{coefficient of overlap - typically in range 0.3 to 1.1}
  \item{\dots}{named arguments passed to \code{\link{Enrm}} and \code{\link{Lambda}} (traps, mask, noccasions, detectfn)}

}

\details{

If sigma is missing and detectfn = `HHN' then sigma is first inferred from the relationship \eqn{\sigma = 100 k \sqrt D} (\eqn{D} in animals per hectare and \eqn{\sigma} in metres). Other detectfn give an error.

A numerical search is then conducted for the value of lambda0 that results in \code{C} expected detections for the given density and design. The calculation takes account of the detector array, the habitat mask and the number of sampling occasions (all specified in the \dots argument - see example).

Only hazard detection functions are supported (`HHN', `HHR', `HEX', `HAN', `HCG'). The default is `HHN'.

}

\value{

A list with one component for each detection parameter.

}

\seealso{

\code{\link{Enrm}}, \code{\link{Lambda}}

}

\examples{

tr <- traps(captdata)
detector(tr) <- "multi"
msk <- make.mask(tr, buffer = 100, type = 'trapbuffer')
getdetectpar(D = 5.48, C = 235, traps = tr, mask = msk, noccasions = 5)

}

\keyword{ datagen }