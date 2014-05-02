##############################################################################
## package 'secrdesign'
## run.scenarios.R
## 2013-02-23, 24, 25, 26, 28
## 2014-02-07, 2014-02-09, 2014-02-10
## currently only for single-session model
## 2014-04-06 drop logfile argument, reorder arguments, add pop.args
## 2014-04-14 fit.models
## 2014-04-26 IHP
## 2014-04-27 automatically wrap mask, pop.args, det.args if not in list form
###############################################################################


###############################################################################
wrapifneeded <- function (args, default) {
    if (any(names(args) %in% names(default)))
        list(args)  ## assume single; wrap
    else
        args
}
###############################################################################
## complete a partial argument list (arg) with default args
## always return a list of arg lists long enough to match max(index)
fullargs <- function (args, default, index) {
    if (is.null(args)) {
        full.args <- list(default)
    }
    else {
        nind <- max(index)
        if (length(args) < nind) stop("too few components in args")
        tmpargs <- vector('list', nind)
        for (i in 1:nind) {
            if (is.character(args[[i]]))
                args[[i]] <- match.arg(args[[i]], default[[i]])
            tmpargs[[i]] <- replace (default, names(args[[i]]), args[[i]])
            if (is.character(tmpargs[[i]]))
                tmpargs[[i]] <- tmpargs[[i]][1]
        }
        full.args <- tmpargs
    }
    full.args
}
###############################################################################

defaultextractfn  <- function(x) {
    if (inherits(x, 'capthist')) {
        ## assume single-session CH
        nmoves <- sum(unlist(sapply(moves(x), function(y) y>0)))
        ## detectors per animal
        dpa <- if (length(dim(x)) == 2)
            mean(apply(abs(x), 1, function(y) length(unique(y[y>0]))))
        else
            mean(apply(apply(abs(x), c(1,3), sum)>0, 1, sum))
        c(n=nrow(x), ndet=sum(abs(x)>0), nmov=nmoves, dpa = dpa)
    }
    else if (inherits(x,'secr'))
        predict(x)
    else
        data.frame()   ## 0 rows, 0 columns
}

###############################################################################
run.scenarios <- function (nrepl,  scenarios, trapset, maskset, xsigma = 4,
    nx = 32, pop.args, det.args, fit = FALSE, fit.args, extractfn = NULL,
    ncores = 1, seed = 123, ...) {

    #--------------------------------------------------------------------------
    onesim <- function (scenario) {
        ## may later allow multi-line scenarios
        if (nrow(scenario) != 1) stop()
        with( scenario, {

            #####################
            ## retrieve data
            grid <- trapset[[trapsindex]]
            mask <- maskset[[maskindex]]
            poparg <- full.pop.args[[popindex]]  ## 2014-04-06
            detarg <- full.det.args[[detindex]]  ## 2014-04-06
            fitarg <- full.fit.args[[fitindex]]

            #####################
            ## override D, core, buffer
            if (poparg$model2D == 'IHP')
                poparg$core <- mask
            else
                poparg$core <- attr(mask, "boundingbox")
            poparg$D <- D*nrepeats
            poparg$buffer <- 0

            #####################
            ## generate population
            pop <- do.call(sim.popn, poparg)

            #####################
            ## form dp for sim.capthist
            ## form par for starting values in secr.fit()
            ## 'par' does not allow for varying link or any non-null model (b, T etc.)
            if (detectfn %in% 14:18) {
                dp <- list(lambda0 = lambda0, sigma = sigma, recapfactor = recapfactor)
                par <- c(log(D), log(lambda0), log(sigma))
            }
            else {
                dp <- list(g0 = g0, sigma = sigma, recapfactor = recapfactor)
                par <- c(log(D), logit(g0), log(sigma))
            }

            #####################
            ## override det args as required
            detarg$traps <- grid
            detarg$popn <- pop
            detarg$detectpar <- dp
            detarg$detectfn <- detectfn
            detarg$noccasions <- noccasions

            #####################
            ## simulate detection
            ## replaced 2014-04-27
            ## CH <- sim.capthist(grid, popn = pop, detectpar = dp,
            ##                   detectfn = detectfn, noccasions = noccasions)
            CH <- do.call(sim.capthist, detarg)
            if (!is.na(nrepeats))
            attr(CH, "n.mash") <- rep(NA,nrepeats)

            #####################
            ## massage results
            if (!fit) {
                extractfn(CH,...)
            }
            else {
                fitarg$capthist <- CH
                fitarg$mask <- mask
                if ((fitarg$start == 'true') | (fitarg$method == 'none')) {
                    ## check to see if simple 'true' values will work
                    ## requires intercept-only model for all parameters
                    model <- eval(fitarg$model)
                    if (!is.list(model)) model <- list(model)
                    vars <- lapply(lapply(model, terms), attr, 'term.labels')
                    if (length(unlist(vars)) == 0) {
                        fitarg$start <- par
                    }
                    else {
                        ## not yet ready for interspersed beta coef
                        warning("using automatic start values")
                        fitarg$start <- NULL
                    }
                }
                fitarg$trace <- FALSE
                fit <- try(do.call(secr.fit, fitarg))
                extractfn(fit,...)
            }
        })
    }
    #--------------------------------------------------------------------------
    runscenario <- function(x) {
        out <- vector('list', nrepl)
        for (r in 1:nrepl) {
            out[[r]] <- onesim(x)
        }
        cat("Completed scenario ", x$scenario, '\n')
        out
    }
    ##--------------------------------------------------------------------------

    ## mainline
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    if (ncores > nrow(scenarios))
        stop ("ncores exceeds number of scenarios")

    ##--------------------------------------------
    ## default extractfn
    if (is.null(extractfn)) {
        extractfn <- defaultextractfn
    }
    ##--------------------------------------------
    ## preprocess inputs
    if (inherits(trapset, 'traps'))   ## otherwise assume already list of traps
        trapset <- list(trapset)
    if (!missing(maskset)) {
        if (inherits(maskset, 'mask'))   ## otherwise assume already list of masks
            maskset <- list(maskset)
        if (is.null(names(maskset)))
            names(maskset) <- paste('mask',1:length(maskset), sep='')
    }
    nk <- length(trapset)
    if (is.null(names(trapset)))
        names(trapset) <- paste('traps',1:nk, sep='')

    dettype <- sapply(trapset, detector)[scenarios$trapsindex]

    OK <- !any((scenarios$nrepeats>1) & (dettype == "single"))
    OK <- if(is.na(OK)) TRUE else OK
    if (!OK)
        warning("single-catch traps violate independence assumption for nrepeats > 1")

    ##---------------------------------------------
    ## allow user changes to default sim.popn arguments
    default.args <- as.list(args(sim.popn))[1:12]
    if (missing(pop.args)) pop.args <- NULL
    pop.args <- wrapifneeded(pop.args, default.args)
    full.pop.args <- fullargs (pop.args, default.args, scenarios$popindex)

    ##---------------------------------------------
    ## allow user changes to default sim.capthist arguments
    default.args <- as.list(args(sim.capthist))[1:13]
    if (missing(det.args)) det.args <- NULL
    det.args <- wrapifneeded(det.args, default.args)
    full.det.args <- fullargs (det.args, default.args, scenarios$detindex)

    ##---------------------------------------------
    ## allow user changes to default secr.fit arguments
    default.args <- as.list(args(secr.fit))[1:21]
    default.args$biasLimit <- NA   ## never check
    default.args$verify <- FALSE   ## never check
    default.args$start <- "true"   ## known values
    if (missing(fit.args)) fit.args <- NULL
    fit.args <- wrapifneeded(fit.args, default.args)
    full.fit.args <- fullargs (fit.args, default.args, scenarios$fitindex)

    ##--------------------------------------------
    ## construct masks as required
    if (missing(maskset)) {
        trapsigma <- scenarios[, c('trapsindex', 'sigma'), drop = FALSE]
        uts <- unique (trapsigma)
        code <- apply(uts, 1, paste, collapse='.')
        scenarios$maskindex <- match(apply(trapsigma,1,paste,collapse='.'), code)
        maskset <- vector('list', nrow(uts))
        for (k in 1:nrow(uts))
            maskset[[k]] <- make.mask(trapset[[uts$trapsindex[k]]], buffer = uts$sigma[k] * xsigma,
                                   type = 'trapbuffer', nx = nx)
    }
    else {
        uts <- NULL
        if (is.null(scenarios$maskindex)) {
            if ((length(maskset) == length(trapset)))
                scenarios[,'maskindex'] <- scenarios$trapsindex
            else if (length(maskset) == 1)
                scenarios[,'maskindex'] <- 1
            else
                stop ("for irregular maskset provide maskindex as a column in scenarios")
        }
        if (max(scenarios$maskindex) > length(maskset))
            stop ("maskindex does not match maskset")
    }

    #--------------------------------------------
    ## override nrepeats and D in scenarios when IHP distribution
    for (i in 1:nrow(scenarios)) {
        pi <- scenarios$popindex[i]
        mi <- scenarios$maskindex[i]
        if (full.pop.args[[pi]]$model2D == 'IHP') {
            avD <- mean (covariates(maskset[[mi]])[,full.pop.args[[pi]]$D])
            scenarios[i, 'nrepeats'] <- 1   ## override
            scenarios[i, 'D'] <- avD
        }
    }

    #--------------------------------------------
    ## run simulations
    tmpscenarios <- split(scenarios, scenarios$scenario)
    if (ncores > 1) {
        ## drop logfile 2014-04-06
        ## clust <- makeCluster(ncores, methods = TRUE, outfile = logfile)
        clust <- makeCluster(ncores, methods = TRUE)
        ## outfile = "" redirects worker output to local console;
        ## does not work on Rgui but works on Rterm and Unix
        clusterSetRNGStream(clust, seed)
        clusterExport(clust,
                      c("onesim", "trapset", "maskset", "fit", "extractfn",
                        "full.pop.args", "full.det.args", "full.fit.args"),
                      environment())
        output <- parLapply(clust, tmpscenarios, runscenario)
        on.exit(stopCluster(clust))
    }
    else {
        set.seed (seed)
        output <- lapply(tmpscenarios, runscenario)
    }
    ##-------------------------------------------
    ## tidy output
    typical <- output[[1]][[1]]
    ## outputtype: secrfit, predicted, coef, numeric
    outputtype <-
        if (inherits(typical, 'secr')) 'secrfit'
        else if (inherits(typical, 'data.frame')) {
                 if (all(c('estimate','SE.estimate','lcl','ucl') %in% names(typical)) &
                     any(c('R.N','E.N') %in% rownames(typical)))
                     'regionN'
                 else if ( all(c('estimate','SE.estimate','lcl','ucl', 'CVn') %in% names(typical)))
                     'derived'
                 else if (all(c('estimate','SE.estimate','lcl','ucl') %in% names(typical)))
                     'predicted'
                 else if (all(c('beta','SE.beta','lcl','ucl') %in% names(typical)))
                     'coef'
                 else 'user'
             }
        else if (inherits(typical, 'capthist'))          ## rawdata
            'capthist'
        else if (is.vector(typical, mode = 'numeric'))   ## usually, summary of unfitted data
            'selectedstatistics'
        else 'user'
    if (outputtype == 'selectedstatistics')
        ## collapse replicates within a scenario into a matrix
        output <- lapply(output, do.call, what = rbind)
    cat("Completed in", round((proc.time() - ptm)[3]/60,3), "minutes \n")
    desc <- packageDescription("secrdesign")  ## for version number
    value <- list (call = cl,
                   version = paste('secrdesign', desc$Version),
                   starttime = starttime,
                   proctime = (proc.time() - ptm)[3],
                   scenarios = scenarios,
                   trapset = trapset,
                   maskset = if (is.null(uts)) maskset else NULL,
                   xsigma = xsigma,
                   nx = nx,
                   pop.args = pop.args,
                   det.args = det.args,
                   fit = fit,
                   fit.args = fit.args,
                   extractfn = extractfn,
                   seed = seed,
                   nrepl = nrepl,
                   output = output,
                   outputtype = outputtype
                   )
    if (outputtype %in% c("secrfit"))
        class(value) <- c("fittedmodels", 'secrdesign', 'list')
    else if (outputtype %in% c("predicted", "derived", "regionN", "coef", "user"))
        class(value) <- c("estimatetables", 'secrdesign', 'list')
    else if (outputtype %in% c("capthist"))
        class(value) <- c("rawdata", 'secrdesign', 'list')
    else if (outputtype %in% c("selectedstatistics"))
        class(value) <- c("selectedstatistics", 'secrdesign', 'list')
    if (outputtype == 'regionN')
        attr(value, 'regionarea') <- sapply(output, function(x) attr(x[[1]], 'regionarea'))
    ## otherwise class remains 'list'

    value
}


########################################################################################

## version of run.scenarios that accepts existing data and
## expands scenarios for multiple model definitions

fit.models <- function (rawdata, fit = FALSE, fit.args, extractfn = NULL,
                           ncores = 1, ...) {

    #--------------------------------------------------------------------------
    onesim <- function (scenario, CH) {
        if (nrow(scenario) != 1) stop()
        with( scenario, {

            ## retrieve data
            grid <- trapset[[trapsindex]]
            mask <- maskset[[maskindex]]
            fitarg <- full.fit.args[[fitindex]]
            if (!fit) {
                extractfn(CH,...)
            }
            else {
                fitarg$capthist <- CH
                fitarg$mask <- mask
                startD <- D*nrepeats

                #####################
                ## form dp for sim.capthist
                ## form par for starting values in secr.fit()
                ## 'par' does not allow for varying link or any non-null model (b, T etc.)
                if (detectfn %in% 14:18) {
                    par <- c(log(startD), log(lambda0), log(sigma))
                }
                else {
                    par <- c(log(startD), logit(g0), log(sigma))
                }

                if ((fitarg$start == 'true') | (fitarg$method == 'none')) {

                    ## check to see if simple 'true' values will work
                    ## requires intercept-only model for all parameters
                    model <- eval(fitarg$model)
                    if (!is.list(model)) model <- list(model)
                    vars <- lapply(lapply(model, terms), attr, 'term.labels')
                    if (length(unlist(vars)) == 0) {
                        fitarg$start <- par
                    }
                    else {
                        ## not yet ready for interspersed beta coef
                        warning("using automatic start values")
                        fitarg$start <- NULL
                    }

                }
                fitarg$trace <- FALSE
                fit <- try(do.call(secr.fit, fitarg))
                extractfn(fit,...)
            }
        })
    }
    #--------------------------------------------------------------------------
    runscenario <- function(x) {
        out <- vector('list', nrepl)
        for (r in 1:nrepl) {
                out[[r]] <- onesim(x, CHlist[[trunc(x$scenario)]][[r]])
        }
        cat("Completed scenario ", x$scenario, '\n')
        out
    }
    ##--------------------------------------------------------------------------

    ## mainline

    if (!inherits(rawdata, "rawdata"))
        stop ("requires rawdata output from run.scenarios()")
    CHlist <- rawdata$output
    nrepl <- rawdata$nrepl
    scenarios <- rawdata$scenarios
    trapset <- rawdata$trapset
    maskset <- rawdata$maskset
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")

    ##--------------------------------------------
    ## default extractfn
    if (is.null(extractfn)) {
        extractfn <- defaultextractfn
    }

    ##---------------------------------------------
    ## allow user changes to default secr.fit arguments
    default.args <- as.list(args(secr.fit))[1:21]
    default.args$biasLimit <- NA   ## never check
    default.args$verify <- FALSE   ## never check
    default.args$start <- "true"   ## known values

    if (missing(fit.args)) fit.args <- NULL
    fit.args <- wrapifneeded(fit.args, default.args)
    nfit <- length(fit.args)
    if (nfit > 1) {
        ## expand scenarios by the required number of different model fits
        scenarios <- scenarios[rep(scenarios$scenario, each = nfit),]
        scenarios$fitindex <- rep(1:nfit, length.out = nrow(scenarios))
        ## assign new unique scenario number by adding decimal fraction
        scenarios$scenario <- scenarios$scenario + scenarios$fitindex /
            10 ^ trunc(log10(nfit)+1)
        rownames(scenarios) <- scenarios$scenario
    }
    full.fit.args <- fullargs (fit.args, default.args, scenarios$fitindex)

    if (ncores > nrow(scenarios))
        stop ("ncores exceeds number of scenarios")

    #--------------------------------------------
    ## run simulations
    tmpscenarios <- split(scenarios, scenarios$scenario)
    if (ncores > 1) {
        ## drop logfile 2014-04-06
        ## clust <- makeCluster(ncores, methods = TRUE, outfile = logfile)
        clust <- makeCluster(ncores, methods = TRUE)
        ## outfile = "" redirects worker output to local console;
        ## does not work on Rgui but works on Rterm and Unix
        clusterExport(clust,
                      c("onesim", "trapset", "maskset", "fit", "extractfn",
                        "full.fit.args", "CHlist"),
                      environment())
        output <- parLapply(clust, tmpscenarios, runscenario)
        on.exit(stopCluster(clust))
    }
    else {
        output <- lapply(tmpscenarios, runscenario)
    }
    ##-------------------------------------------
    ## tidy output
    typical <- output[[1]][[1]]
    ## outputtype: secrfit, predicted, coef, numeric
    outputtype <-
        if (inherits(typical, 'secr')) 'secrfit'
        else if (inherits(typical, 'data.frame')) {
                 if (all(c('estimate','SE.estimate','lcl','ucl') %in% names(typical)) &
                     any(c('R.N','E.N') %in% rownames(typical)))
                     'regionN'
                 else if ( all(c('estimate','SE.estimate','lcl','ucl', 'CVn') %in% names(typical)))
                     'derived'
                 else if (all(c('estimate','SE.estimate','lcl','ucl') %in% names(typical)))
                     'predicted'
                 else if (all(c('beta','SE.beta','lcl','ucl') %in% names(typical)))
                     'coef'
                 else 'user'
             }
        else if (inherits(typical, 'capthist'))
            'capthist'
        else if (is.vector(typical, mode = 'numeric'))
            'selectedstatistics'
        else 'user'
    if (outputtype == 'selectedstatistics')
        ## collapse replicates within a scenario into a matrix
        output <- lapply(output, do.call, what = rbind)
    cat("Completed in", round((proc.time() - ptm)[3]/60,3), "minutes \n")
    desc <- packageDescription("secrdesign")  ## for version number
    value <- list (call = cl,
                   version = paste('secrdesign', desc$Version),
                   starttime = starttime,
                   proctime = (proc.time() - ptm)[3],
                   scenarios = scenarios,
                   trapset = trapset,
                   maskset = maskset,
                   xsigma = rawdata$xsigma,
                   nx = rawdata$nx,
                   pop.args = rawdata$pop.args,
                   det.args = rawdata$det.args,
                   fit = fit,
                   fit.args = fit.args,
                   extractfn = extractfn,
                   seed = rawdata$seed,
                   nrepl = rawdata$nrepl,
                   output = output,
                   outputtype = outputtype
                   )
    if (outputtype %in% c("secrfit"))
        class(value) <- c("fittedmodels", 'secrdesign', 'list')
    else if (outputtype %in% c("predicted", "derived", "regionN", "coef", "user"))
        class(value) <- c("estimatetables", 'secrdesign', 'list')
    else if (outputtype %in% c("capthist"))
        class(value) <- c("rawdata", 'secrdesign', 'list')
    else if (outputtype %in% c("selectedstatistics"))
        class(value) <- c("selectedstatistics", 'secrdesign', 'list')
    if (outputtype == 'regionN')
        attr(value, 'regionarea') <- sapply(output, function(x) attr(x[[1]], 'regionarea'))
    ## otherwise class remains 'list'

    value
}
