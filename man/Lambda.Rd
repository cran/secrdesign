\name{Lambda}
\alias{Lambda}
%\alias{Lambdak}
\alias{Enrm}
\alias{minnrRSE}
\alias{Qpm}
\alias{En2}

\title{Expected Detections}

\description{

Compute the expected number of detections as a function of location (\code{Lambda}), and the expected total numbers of individuals \eqn{n}, recaptures \eqn{r} and movements \eqn{m} for a population sampled with an array of detectors (\code{Enrm}) or the number of individuals detected at two or more detectors (\code{En2}). 

}

\usage{

Lambda(traps, mask, detectpar, noccasions, detectfn = c("HHN", "HHR", "HEX", 
    "HAN", "HCG", 'HN', 'HR', 'EX'))
% Lambdak(D, traps, mask, detectpar, detectfn = c("HHN", "HHR", "HEX", 
%    "HAN", "HCG", 'HN', 'HR', 'EX'))
Enrm(D, ...)

minnrRSE(D, ..., CF = 1.0, distribution = c("poisson","binomial"))

En2(D, traps, mask, detectpar, noccasions, detectfn = c("HHN", "HHR", "HEX", 
    "HAN", "HCG", "HN", "HR", "EX"))

Qpm(D, traps, mask, detectpar, noccasions, detectfn = c("HHN", "HHR", "HEX", 
    "HAN", "HCG", "HN", "HR", "EX"))

}

\arguments{
  \item{traps}{ \code{\link[secr]{traps}} object }
  \item{mask}{ \code{\link[secr]{mask}} object }
  \item{detectpar}{ a named list giving a value for each parameter of detection function }
  \item{noccasions}{integer number of sampling occasions}
  \item{detectfn}{integer code or character string for shape of detection function 
      -- see \link[secr]{detectfn}}
  \item{D}{population density animals / hectare; may be scalar or vector of length \code{nrow(mask)}}
  \item{\dots}{arguments passed to \code{Lambda}}
  \item{CF}{numeric correction factor}
  \item{distribution}{character distribution of \eqn{n}}
}

\details{

The detector attribute of \code{traps} may be `multi', `proximity' or `count'. It is assumed that detectpar and detector type do not differ among occasions.

The calculation is based on an additive hazard model. If \code{detectfn} is not a hazard function (`HHN', `HEX', `HHR', `HAN' and `HCG') then an attempt is made to approximate one of the hazard functions (HN -> HHN, HR -> HHR, EX -> HEX). The default is `HHN'.

For hazard function \eqn{\lambda(d)} and \eqn{S} occasions, we define \eqn{\Lambda(x) = \sum_s \sum_k \lambda(d_k(x))}{\Lambda = \sum over s \sum over k \lambda(d_k(x)) }.

Formulae for expected counts are given in \href{https://www.otago.ac.nz/density/pdfs/secrdesign-Enrm.pdf}{secrdesign-Enrm.pdf}.

\code{minnrRSE} has mostly the same inputs as \code{Enrm} but returns sqrt(CF/min(n,r)). 
The correction factor CF may be used to adjust for systematic bias (e.g., for a line of
detectors CF = 1.4 may be appropriate). The default \code{distribution = 'poisson'} is 
for Poisson-distributed \eqn{N} and \eqn{n}. To adjust the prediction for fixed 
\eqn{N} (binomial \eqn{n}) use \code{distribution = 'binomial'} 
(see \href{https://www.otago.ac.nz/density/pdfs/secrdesign-tools.pdf}{secrdesign-tools.pdf} Appendix 2).

From 2.7.0, the first argument of \code{minnrRSE} may also be the output from \code{\link{GAoptim}}.

\code{En2} is defined for detectors `multi', `proximity' and `count'.

\code{Qpm} returns the optimisation criteria \eqn{Q_p} and \eqn{Q_{p_m}} of Dupont et al. (2021), defined only for `proximity' and `count' detectors. The criteria are mask-dependent, and \code{En2} is generally preferred. For `proximity' and `count' detectors the following expressions give the same result:

\code{En2(D, trp, msk, dp)}

\code{Qpm(D, trp, msk, dp) * maskarea(msk) * D}

given constant density `D', detectors `trp', mask `msk' and detection parameters `dp'.

}

\value{

\code{Lambda} --
    \link[secr]{mask} object with covariates `Lambda' (\eqn{\Lambda(x)}{\Lambda(x)}), `sumpk' and `sumq2' (intermediate values for computation of expected counts)
%\code{Lambdak} --
%    Vector of with single-occasion hazard at each detector integrated over \link{mask} object.

\code{Enrm} --
    numeric vector of length 3, the values of E(\eqn{n}), E(\eqn{r}) and E(\eqn{m})

\code{minnrRSE} -- rule-of-thumb RSE(D-hat) Efford and Boulanger (2019)

\code{En2} -- numeric vector comprising the values E(\eqn{n}) and E(number of animals detected at 2 or more sites)

\code{Qpm} -- numeric vector comprising the criteria \eqn{Q_p} and \eqn{Q_{p_m}} of Dupont et al. (2021)


}

\references{

Dupont, G., Royle, J. A., Nawaz, M. A. and Sutherland, C. (2021) Optimal sampling
design for spatial capture--recapture. \emph{Ecology} \bold{102} e03262.

Efford, M. G., and Boulanger, J. (2019) Fast evaluation of study designs for 
spatially explicit capture--recapture. 
\emph{Methods in Ecology and Evolution}, \bold{10}, 1529--1535. \doi{10.1111/2041-210X.13239}

}

\seealso{

\code{\link{getdetectpar}}, 
\code{\link{optimalSpacing}}, 
\code{\link{scenarioSummary}}, 
\code{\link{GAoptim}}

}

\examples{

tr <- traps(captdata)
detector(tr) <- "multi"
msk <- make.mask(tr, buffer = 100, type = 'trapbuffer')

L <- Lambda(tr, msk, list(lambda0 = 0.2, sigma = 20), 5)
nrm <- Enrm(D = 5, tr, msk, list(lambda0 = 0.2, sigma = 20), 5)
nrm

En2(D = 5, tr, msk, list(lambda0 = 0.2, sigma = 20), 5)

plot(L, cov = "Lambda", dots = FALSE)
plot(tr, add = TRUE)
mtext(side = 3,  paste(paste(names(nrm), round(nrm,1)), collapse = ", "))

}

\keyword{ manip }