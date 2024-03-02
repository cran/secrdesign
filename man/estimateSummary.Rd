\name{estimateSummary}
\alias{estimateArray}
\alias{estimateSummary}

\title{Direct summary of estimate tables}

\description{

An alternative approach to summarising output from \code{run.scenarios}
(cf \code{\link{summary.estimatetables}}). \code{estimateSummary} is 
especially useful when \code{extractfn = predict} or \code{extractfn = coef},
and all scenarios have group structure.

}

\usage{

estimateArray(object)

estimateSummary(object, parameter = "D", statistics = c("true", "nvalid", 
    "EST", "seEST", "RB", "seRB", "RSE", "RMSE", "rRMSE", "COV"), true,
    validrange = c(0, Inf), checkfields = c('estimate','SE.estimate'),
    format = c('data.frame', 'list'), cols = NULL)

}

\arguments{

  \item{object}{secrdesign object of class "estimatetables"}
  \item{parameter}{character name of parameter (row in estimate table)}
  \item{statistics}{character choice of outputs}
  \item{true}{numeric vector of true values, one per scenario and group}
  \item{validrange}{numeric allowed for estimates or other checkfields}
  \item{checkfields}{character choice of columns in each estimate table that 
  will be checked against validrange}
  \item{format}{character choice of output}
  \item{cols}{indices of scenario columns to include when format = "data.frame"}
  
}

\details{

When 'predict(fittedmodel)' in \code{\link{run.scenarios}} generates more than 
one estimate table (i.e. when the model uses groups, mixture classes or 
multiple sessions), the default extract function retains only the first. 
This is often OK, but it can be frustrating if you care about 
group- or session-specific estimates.

The alternative is to use 'predict' as the run.scenarios extractfn, which 
retains all estimate tables. This requires a different function for 
summarisation; \code{estimateSummary} will suffice for many purposes.

\code{estimateSummary} internally calls \code{estimateArray} to pre-process 
the output from run.scenarios.

The code should be examined for the precise definition of each statistic.

True parameter values are required for RB, RMSE and COV, and these are computed 
even if later dropped from the output. If provided, the argument \code{true} 
should have length equal to the number of parameter tables in each replicate, 
i.e. (number of scenarios) * (number of groups), ordered by scenario. 
Otherwise, true values will be taken from rows of the data frame object$scenarios.

Replicates are rejected (set to NA) if any \code{checkfields} 
falls outside \code{validrange}.

Output statistics `EST', `RB', and `RSE' are means 
across replicates, and `seEST', `seRB' the corresponding standard errors.

The output list may optionally be formatted as a data.frame with pre-pended 
columns from object$scenarios. Set \code{cols} to 0 or NULL for no scenario 
columns.

\code{cols} defaults to \code{c("scenario", "group")} if groups are present and \code{"scenario"} otherwise. 

}

\value{

\code{estimateArray} --
array with dimensions (Parameter, statistic, Group, Scenario, Replicate)

\code{estimateSummary} --

If groups present and format = "list" - a list of matrices (group x scenario), one for each statistic:

\item{true.X}{true value of parameter (X)}
\item{nvalid}{number of valid replicates used in later summaries}
\item{EST}{mean of parameter estimates}
\item{seEST}{standard error of estimates (across replicates)}
\item{RB}{relative bias}
\item{seRB}{standard error of replicate-specific RB (across replicates)}
\item{RSE}{relative standard error (SE.estimate/estimate)}
\item{RMSE}{root mean squared error}
\item{rRMSE}{RMSE/true.X}
\item{COV}{coverage of confidence intervals (usually 95\% intervals).}

If groups absent and format = "list" - a list of vectors (one element per scenario) with statistics as above.

If format = "data.frame" - a data frame with rows corresponding to group x scenario (or session x scenario) combinations and columns corresponding to statistics as above.

}

\note{

These functions were introduced in version 2.8.1. They may change in later versions. 
The default format was changed to `data.frame' in 2.8.3.

Results may be confusing when scenarios have group structure and groups 
are not used in the fitted model.

It is assumed that all scenarios (and all models in a multi-model fit) report the same parameters.

}

\seealso{

\code{\link{run.scenarios}},
\code{\link{header}},
\code{\link{summary.estimatetables}}

}

\examples{

\donttest{

# 2-scenario, 2-group simulation
scen8 <- make.scenarios (D = 8, g0 = 0.3, sigma = 30, 
    noccasions = c(4,8), groups = c('F','M'))
    
# replace density and sigma values of males to make it interesting 
male <- scen8$group == 'M'
scen8$D[male] <- 4
scen8$sigma[male] <- 40

grid <- make.grid(8, 8, spacing = 30)
mask <- make.mask(grid, buffer = 160, type = 'trapbuffer')

old <- options(digits = 3)
setNumThreads(2)

#--------------------------------------------------------------------------
# run a few simulations

# model groups
sims <- run.scenarios(10, scen8, trapset = grid, fit = TRUE, 
    fit.args = list(model = list(D~g, g0~1, sigma~g), groups = 'group'),
    extractfn = predict, maskset = mask)

# format as list, selecting statistics
# default summary uses true = c(8,4,8,4)  
estimateSummary(sims, 'D', c("true", "nvalid", "EST", "RB", "seRB"))

# format as data.frame by scenario and group, all statistics
estimateSummary(sims, 'D',  format = 'data.frame')

#--------------------------------------------------------------------------
# try with default extractfn (single table per replicate, despite groups)
sims2 <- run.scenarios(10, scen8, trapset = grid, fit = TRUE, 
     maskset = mask)

# Fails with "Error in estimateSummary(sims2, "D") : incongruent 'true'""
# estimateSummary(sims2, 'D')

# OK if manually provide scenario-specific true density
estimateSummary(sims2, 'D', true = c(12,12))

# reformat by scenario
estimateSummary(sims2, 'D', true = c(12,12), format = 'data.frame')

# compare standard summary
summary(sims2)$OUTPUT

#--------------------------------------------------------------------------

# multiple estimate tables also arise from multi-session simulations
# argument 'true' must be specified manually
# interpret with care: sessions are (probably) not independent
# this example uses the previous grid and mask

scen9 <- make.scenarios (D = 8, g0 = 0.3, sigma = 30, noccasions = 5)
poparg <- list(nsessions = 3, details = list(lambda = 1.2))  # for sim.popn
detarg <- list(renumber = FALSE)                             # for sim.capthist
fitarg <- list(model = D~Session)                            # for secr.fit

sims3 <- run.scenarios(5, scen9, trapset = grid, fit = TRUE, 
    maskset = mask, pop.args = poparg, det.args = detarg,
    fit.args = fitarg, extractfn = predict)
    
estimateSummary(sims3, parameter = 'D', format = 'data.frame', 
    true = 8 * 1.2^(0:2))
#--------------------------------------------------------------------------

# extractfn = coef results in a single estimate table per replicate,
# so the usual summary method is sufficent. For completeness we show 
# that estimateSummary can also be used. Coefficients are often negative,
# so relative values (e.g., RB, RSE) may be meaningless.

sims4 <- run.scenarios(5, scen9, trapset = grid, fit = TRUE, 
    maskset = mask, pop.args = poparg, det.args = detarg,
    fit.args = fitarg, extractfn = coef)
    
estimateSummary(sims4, parameter = 'D', c("nvalid", "EST", "seEST", "RMSE", "COV"), 
    format = 'data.frame', true = log(8), checkfields = 'beta', 
    validrange = log(c(2,20)))

estimateSummary(sims4, parameter = 'D.Session', c("nvalid", "EST", "seEST", 
    "RMSE", "COV"), format = "data.frame", true = log(1.2), checkfields = "beta", 
    validrange = log(c(0.5,2)))
#--------------------------------------------------------------------------

options(old)

}

}
