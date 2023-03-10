\name{run.scenarios}
\alias{run.scenarios}
\alias{fit.models}

\title{Simulate Sampling Designs}

\description{

This function performs simulations to predict the precision of density and other
estimates from simple 1-session SECR designs. Scenarios are specified
via an input dataframe that will usually be constructed with
\code{\link{make.scenarios}}. Each scenario comprises an index to a detector layout,
the number of sampling occasions, and specified density (D) and detection
parameters (usually \eqn{g_0} and \eqn{\sigma}).

Detector layouts are provided in a separate list \code{trapset}. This
may comprise an actual field design input with \code{\link{read.traps}} or
`traps' objects constructed with \code{\link{make.grid}} etc., as in the
Examples. Even a single layout must be presented as a component of a
list (e.g., \code{list(make.grid())}).

Alternative approaches are offered for predicting precision. Both start
by generating a pseudorandom dataset under the design using the
parameter values for a particular scenario. The first estimates the
parameter values and their standard errors from each dataset by
maximizing the full likelihood, as usual in \code{secr.fit}. The second
takes the short cut of computing variances and SE from the Hessian
estimated numerically at the known expected values of the parameters,
without maximizing the likelihood. Set \code{method = "none"} in fit.args 
for this shortcut.

}

\usage{

run.scenarios(nrepl, scenarios, trapset, maskset, xsigma = 4, nx = 32,
    pop.args, CH.function = c("sim.capthist", "simCH"), det.args, 
    fit = FALSE, fit.function = c("secr.fit", "ipsecr.fit"), 
    fit.args, chatnsim, extractfn = NULL, multisession = FALSE,
    ncores = NULL, byscenario = FALSE, seed = 123, trap.args, ...)

fit.models(rawdata, fit = FALSE, fit.function = c("secr.fit", "ipsecr.fit"), 
    fit.args, chatnsim, extractfn = NULL, ncores = NULL, byscenario = FALSE,
    scen, repl, ...)
}

\arguments{
  \item{nrepl}{integer number of replicate simulations}
  \item{scenarios}{dataframe of simulation scenarios}
  \item{trapset}{secr traps object or a list of traps objects or functions}
  \item{maskset}{secr mask object or a list of mask objects (optional)}
  \item{xsigma}{numeric buffer width as multiple of sigma (alternative
  to maskset)}
  \item{nx}{integer number of cells in mask in x direction  (alternative
  to maskset)}
  \item{pop.args}{list of named arguments to
  \code{\link[secr]{sim.popn}} (optional)} 
  \item{CH.function}{character name of function to simulate capthist}
  \item{det.args}{list of named arguments to
  \code{\link[secr]{sim.capthist}} (optional)} 
  \item{fit}{logical; if TRUE a model is fitted with \code{fit.function}, otherwise
  data are generated but no model is fitted \cr (see also Design-only statistics in Details)}
  \item{fit.function}{character name of function to use for model fitting}
  \item{fit.args}{list of named arguments to \code{\link[secr]{secr.fit}} or \code{\link[ipsecr]{ipsecr.fit}} (optional)}
  \item{chatnsim}{integer number of simulations for overdispersion of mark-resight models}
  \item{extractfn}{function to extract a vector of statistics from secr model}
  \item{multisession}{logical; if TRUE groups are treated as additional sessions}
  \item{ncores}{integer number of cores for parallel processing or NULL}
  \item{byscenario}{logical; if TRUE then each scenario is sent to a different core}
  \item{seed}{integer pseudorandom number seed}
  \item{trap.args}{list of arguments for trapset components if using function option }
  \item{\dots}{other arguments passed to extractfn}
  \item{rawdata}{`rawdata' object from previous call to \code{run.scenarios}}
  \item{scen}{integer vector of scenario subscripts}
  \item{repl}{integer vector of subscripts in range 1:nrepl}
  }

\details{

  Designs are constructed from the trap layouts in \code{trapset}, the
  numbers of grids in \code{ngrid}, and the numbers of sampling
  occasions (secondary sessions) in \code{noccasions}. These are
  \emph{not} crossed: the number of designs is the maximum length of any
  of these arguments. Any of these arguments whose length is less than
  the maximum will be replicated to match.

  \code{pop.args} is used to customize the simulated population
  distribution. It will usually comprise a single list, but may be a
  list of lists (one per popindex value in scenarios).

  \code{det.args} may be used to customize some aspects of the detection
  modelling in \code{sim.capthist}, but not \code{traps, popn, detectpar,
  detectfn}, and \code{noccasions}, which are controlled directly by the
  scenarios. It will usually comprise a single list, but may be a list
  of lists (one per detindex value in scenarios).

  \code{fit.args} is used to customize the fitted model; it will usually
  comprise a single list. If you are interested in precision alone, use
  \code{fit.args=list(method = 'none')} to obtain variance estimates
  from the hessian evaluated at the parameter estimates. This is much
  faster than a complete model fit, and usually accurate enough.
  
  If no \code{extractfn} is supplied then a default is used - see
  Examples. Replacement functions should follow this pattern i.e. test
  for whether the single argument is an secr object, and if not supply a
  named vector of NA values of the correct length.
  
  Using \code{extractfn = summary} has the advantage of allowing both model fits and raw statistics to be extracted from one set of simulations. However, this approach requires an additional step to retrieve the desired numeric results from each replicate (see \code{\link{count.summary}} and \code{\link{predict.summary}}).

\subsection{Parallel processing}{  
  If \code{byscenario = TRUE} then by default each scenario will be run in a separate worker
process using \code{parLapply} from \pkg{parallel} (see also \link[secr]{Parallel}). The number of scenarios should not exceed the available number of cores (set by the 'ncores' argument or a prior call to `setNumThreads`).

If \code{byscenario = FALSE} then from \pkg{secrdesign} 2.6.0 onwards the usual multithreading of \pkg{secr} 4.5 is applied. The number of cores should usually be preset with `setNumThreads`. If \code{ncores} is provided then the environment variable RCPP_PARALLEL_NUM_THREADS is reset. The default behaviour of the fitting functions (secr.fit, ipsecr.fit) is to use this value (unless specified in fit.args).

  When `byscenario = TRUE` the L'Ecuyer pseudorandom generator is used with a separate random
  number stream for each core (see \code{\link{clusterSetRNGStream}}).

For \code{ncores > 1} it pays to keep an eye on the processes from the
  Performance page of Windows Task Manager (<ctrl><alt><del>), or `top' in
  linux OS. If you interrupt \code{run.scenarios} (<Esc> from Windows)
  you may occasionally find some processes do not terminate and have to
  be manually terminated from the Task Manager - they appear as
  Rscript.exe on the Processes page.
  
}

\subsection{Alternate functions for simulation and fitting}{

The default is to use functions \code{\link{sim.capthist}} and \code{\link{secr.fit}} from \pkg{secr}. Either may be substituted by the corresponding function (\code{\link[ipsecr]{simCH}} or \code{\link[ipsecr]{ipsecr.fit}}) from package \pkg{ipsecr} if that has been installed.

}

\subsection{Design-only statistics}{

Designs for distance sampling were evaluated by Fewster and Buckland (2004) by computing statistics from simulated detections without fitting a model to estimate the detection parameters. An analogous procedure for SECR is implemented by setting \code{fit = 'design'}. A new default extractfn (designextractfn) computes the effective sampling area with the \pkg{secr} function \code{\link{pdot}} and returns a vector of values -

\tabular{lll}{
  n   \tab \tab number of individuals detected \cr
  r   \tab \tab number of recaptures \cr
  esa \tab \tab effective sampling area, given the known detection parameters \cr
  D   \tab \tab D = n/esa \cr
}

The resulting simulation object is of type 'selectedstatistics' for which the summary method works as usual.

A similar effect may be achieved by providing a custom extractfn and passing arguments to it via the dots argument of \code{run.scenarios}.

}

\subsection{Miscellaneous}{

From 2.2.0, two or more rows in \code{scenarios} may share the same scenario number. This is used to generate multiple population subclasses (e.g. sexes) differing in density and/or detection parameters. If \code{multisession = TRUE} the subclasses become separate sessions in a multi-session capthist object (this may require a custom \code{extractfn}). \code{multisession} is ignored with a warning if each scenario row has a unique number.
  
From 2.7.0, each component of `trapset' may be a function that constructs a detector layout. This allows layouts to be constructed dynamically at the time each capthist is generated; arguments of each function are provided in the `trap.args' list. The primary purpose is to allow systematic grids, laceworks etc. to be constructed with a unique random origin for each replicate. The `maskset' argument must be provided - it should cover all potential layouts, regardless of origins.
  
In \code{fit.models} the arguments \code{scen} and \code{repl} may be used to select a subset of datasets for model fitting.
  
}
  
\subsection{Mark-resight}{  \code{chatnsim} controls an additional quasi-likelihood model step to adjust for overdispersion of sighting counts. No adjustment happens when \code{chatnsim = 0}; otherwise \code{abs(chatnsim)} gives the number of simulations to perform to estimate overdispersion. If \code{chatnsim < 0} then the quasilikelihood is used only to re-estimate the variance at the previous MLE (method = "none").
}


 \subsection{Further processing}{A summary method is provided (see
  \code{\link{summary.secrdesign}}). It is usually necessary to process
  the simulation results further with \code{\link{predict.fittedmodels}}
  and/or \code{\link{select.stats}} before summarization.
  }
  
}

\value{

  An object of class (x, `secrdesign', `list'), where x is one of
  `fittedmodels', `estimatetables', `selectedstatistics' or `rawdata',
  with components

  \item{call}{function call}
  \item{version}{character string including the software version number}
  \item{starttime}{character string for date and time of run}
  \item{proctime}{processor time for simulations, in seconds}
  \item{scenarios}{dataframe as input}
  \item{trapset}{list of trap layouts as input}
  \item{maskset}{list of habitat masks (input or generated)}
  \item{xsigma}{from input}
  \item{nx}{from input}
  \item{pop.args}{from input}
  \item{CH.function}{from input}
  \item{det.args}{from input}
  \item{fit}{from input}
  \item{fit.function}{from input}
  \item{fit.args}{from input}
  \item{extractfn}{function used to extract statistics from each
    simulation}
  \item{seed}{from input}
  \item{nrepl}{from input}
  \item{output}{list with one component per scenario }
  \item{outputtype}{character code - see vignette}

  If \code{fit = FALSE} and \code{extractfn = identity} the result is of
  class (`rawdata', `secrdesign', `list'). This may be used as input to
  \code{fit.models}, which interprets each model specification in
  \code{fit.args} as a new `sub-scenario' of each input scenario
  (i.e. all models are fitted to every dataset). The output
  possibilities are the same as for \code{run.scenarios}.
  
 If subclasses have been defined (i.e. scenarios has multiple rows with the same scenario ID), each simulated capthist object has covariates with a character-valued column named "group" ("1", "2" etc.) (there is also a column "sex" generated automatically by \code{sim.popn}).
  
}

\note{

100 ha = 1 km^2.

fit.function = 'openCR.fit' was deprecated from 2.5.8 and has been removed.
  
}

\author{
Murray Efford
}


\references{

  Fewster, R. M. and Buckland, S. T. 2004. Assessment of distance
  sampling estimators. In: S. T. Buckland, D. R. Anderson,
  K. P. Burnham, J. L. Laake, D. L. Borchers and L. Thomas (eds)
  \emph{Advanced distance sampling}. Oxford University Press, Oxford,
  U. K. Pp. 281--306.
  
}


\seealso{
  \code{\link{select.stats}},
  
  \code{\link{summary.secrdesign}},
  
  \code{\link{summary.estimatetables}},
  
  \code{\link{summary.selectedstatistics}},
  
  \code{\link{estimateSummary}}
  

  Miscellaneous --
  
  \code{\link{predict.fittedmodels}},
  
  \code{\link{scenarioSummary}}, 
  
  \code{\link{count.summary}},
  
  \code{\link{predict.summary}}
  
  
  \pkg{secr} functions used internally --
  
  \code{\link[secr]{sim.popn}},
  
  \code{\link[secr]{sim.capthist}},
  
  \code{\link[secr]{secr.fit}}
  
  
  To combine output --
  
  \code{\link{rbind.estimatetables}},
  
  \code{\link{rbind.selectedstatistics}},
  
  \code{\link{c.estimatetables}},
  
  \code{\link{c.selectedstatistics}}
  
}

\examples{

## Simple example: generate and summarise trapping data
## at two densities and for two levels of sampling frequency
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2, noccasions =
    c(5,10))
traps1 <- make.grid()   ## default 6 x 6 trap grid
tmp1 <- run.scenarios(nrepl = 20, trapset = traps1, scenarios = scen1,
    fit = FALSE)
summary(tmp1)

\dontrun{

setNumThreads(7)

##########################################
# new summary method (secrdesign >= 2.8.1)
# assumes fit = TRUE, extractfn = predict

tmp2 <- run.scenarios(nrepl = 10, trapset = traps1, scenarios = scen1,
    fit = TRUE, extractfn = predict)
estimateSummary(tmp2, format = "data.frame", 
    cols = c('scenario', 'noccasions'))

###########################
## 2-phase example
## first make and save rawdata
scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
traps1 <- make.grid()   ## default 6 x 6 trap grid
tmp1 <- run.scenarios(nrepl = 20, trapset = traps1, scenarios = scen1,
    fit = FALSE, extractfn = identity)

## review rawdata
summary(tmp1)

## then fit and summarise models
tmp2 <- fit.models(tmp1, fit.args = list(list(model = g0~1),
    list(model = g0~T)), fit = TRUE)
summary(tmp2)
###########################

## Construct a list of detector arrays
## Each is a set of 5 parallel lines with variable between-line spacing;
## the argument that we want to vary (spacey) follows nx, ny and spacex
## in the argument list of make.grid().

spacey <- seq(2000,5000,500)
names(spacey) <- paste('line', spacey, sep = '.')
trapset <- lapply(spacey, make.grid, nx = 101, ny = 5, spacex = 1000,
    detector = 'proximity')

## Make corresponding set of masks with constant spacing (1 km)
maskset <- lapply(trapset, make.mask, buffer = 8000, spacing = 1000,
    type = 'trapbuffer')

## Generate scenarios
scen <- make.scenarios (trapsindex = 1:length(spacey), nrepeats = 8,
    noccasions = 2, D = 0.0002, g0 = c(0.05, 0.1), sigma = 1600, cross = TRUE)

## RSE without fitting model
sim <- run.scenarios (50, scenarios = scen, trapset = trapset, maskset = maskset,
    fit = TRUE, fit.args = list(method = 'none'), seed = 123)

## Extract statistics for predicted density
sim <- select.stats(sim, parameter = 'D')

## Plot to compare line spacing
summ <- summary (sim, type='array',  fields = c('mean','lcl','ucl'))$OUTPUT
plot(0,0,type='n', xlim=c(1.500,5.500), ylim = c(0,0.36), yaxs = 'i',
    xaxs = 'i', xlab = 'Line spacing  km', ylab = 'RSE (D)')
xv <- seq(2,5,0.5)
points(xv, summ$mean[,1,'RSE'], type='b', pch=1)
points(xv, summ$mean[,2,'RSE'], type='b', pch=16)
segments(xv, summ$lcl[,1,'RSE'], xv, summ$ucl[,1,'RSE'])
segments(xv, summ$lcl[,2,'RSE'], xv, summ$ucl[,2,'RSE'])
legend(4,0.345, pch=c(1,16), title = 'Baseline detection',
    legend = c('g0 = 0.05', 'g0 = 0.1'))
}

}

\keyword{ Datagen }

%## the 20:20 rule
% recapt <- apply(summary(out8nsigma,c('n','ncapt'))$mean, 1:2, diff)
% plot(recapt, meanRSE, xlim=c(0,90), ylim=c(0,0.7),
%     xaxs = 'i', yaxs = 'i',
%     xlab = 'Number of recaptures', ylab = 'RSE(D-hat)')
% lines(c(0,20,20),c(0.2,0.2,0), lty=2)
% lines(x <- 0:85, 1/x^0.5)
% mtext (side=3, line = 1, 'The 20:20 rule and 1/nrecapt^0.5', cex=0.8)
