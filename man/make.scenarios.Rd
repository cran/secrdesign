\name{make.scenarios}
\alias{make.scenarios}

\title{Construct Scenario Data Frame}

\description{

  This function prepares a dataframe in which each row specifies a
  simulation scenario. The dataframe is used as input to \code{\link{run.scenarios}}.

}

\usage{

make.scenarios(trapsindex = 1, noccasions = 3,  nrepeats = 1, D, g0, sigma, lambda0,
detectfn = 0, recapfactor = 1, popindex = 1, detindex = 1, fitindex = 1, groups,
crosstraps = TRUE)

}

\arguments{ 
  \item{trapsindex}{integer vector determining the traps object to use}
  
  \item{noccasions}{integer vector for the number of sampling occasions}
  
  \item{nrepeats}{integer vector of multipliers for D (see Details)}
  
  \item{D}{numeric vector of values for the density parameter (animals /
    hectare)}
  
  \item{g0}{numeric vector of values for the g0 parameter}
  
  \item{sigma}{numeric vector of values for the sigma parameter (m)}
  
  \item{lambda0}{numeric vector of values for the lambda0 parameter}
  
  \item{detectfn}{vector of valid detection function codes (numeric or
    character)}
  
  \item{recapfactor}{numeric vector of values for recapfactor
    (\code{\link{sim.capthist}})}
  
  \item{popindex}{integer vector determining which population model is used}
  
  \item{detindex}{integer vector determining which detection options are used}
  
  \item{fitindex}{integer vector determining which model is fitted}
  
  \item{groups}{character vector of group labels (optional)}
  
  \item{crosstraps}{logical; if TRUE the output includes all
    combinations of   
  \code{trapsindex}, \code{noccasions} and \code{nrepeats}}
}

\details{

  The index in \code{trapsindex} is used in \code{\link{run.scenarios}} to
  select particular detector arrays from the list of arrays provided as
  an argument to that function.
  
  The function generates all combinations of the given parameter values
  using \code{\link[base]{expand.grid}}. By default, it also generates
  all combinations of the parameters with \code{trapsindex} and the
  number of sampling occasions. If \code{crosstraps} is FALSE then
  \code{trapsindex}, \code{noccasions}, and \code{nrepeats} are merely
  used to fill in these columns in the output dataframe.

  The argument \code{lambda0} replaces \code{g0} for the hazard detection
  functions 14--18 (\code{\link[secr]{detectfn}}). 
  
  Designs may use multiple detector arrays with the same internal
  geometry (e.g., number and spacing of traps). The number of such
  arrays is varied with the \code{nrepeats} argument. For example, you
  may compare designs with many small arrays or a few large ones. In
  practice, \code{run.scenarios} simulates a single layout with density 
  D * nrepeats. This shortcut is not appropriate when
  animals compete for traps (detector = `single').

  \code{fitindex} allows a choice of different models when the argument
  \code{fit.args} of \code{\link{run.scenarios}} is a compound list.
  
  If \code{groups} is provided each scenario is replicated to the length of \code{groups} and a column `group' is added.
  
}

\value{

  Dataframe with one row per scenario (or sub-scenario) and the columns
  \item{scenario}{a number identifying the scenario}
  \item{group}{(optional)}
  \item{trapsindex}{}
  \item{noccasions}{}
  \item{nrepeats}{}
  \item{D}{}
  \item{g0}{or lambda0}
  \item{sigma}{}
  \item{detectfn}{see \code{\link[secr]{detectfn}}; always numeric}
  \item{recapfactor}{}
  \item{popindex}{}
  \item{detindex}{}
  \item{fitindex}{}
  
  An attribute `inputs' is saved for possible use in \code{\link{make.array}}.
  
}

\seealso{
  
  \code{\link{run.scenarios}}, \code{\link{scenarioSummary}}, \code{\link[secr]{sim.capthist}}
  
}

\examples{

make.scenarios(trapsindex = 1, nrepeats = 1, D = c(5,10), sigma = 25,
g0 = 0.2)

}

\keyword{ manip }
