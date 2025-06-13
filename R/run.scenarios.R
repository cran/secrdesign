#############################################################################
## package 'secrdesign'
## run.scenarios.R
## 2013-02-23, 24, 25, 26, 28
## 2014-02-07, 2014-02-09, 2014-02-10, 2017-05-06
## currently only for single-session model
## 2014-04-06 drop logfile argument, reorder arguments, add pop.args
## 2014-04-14 fit.models
## 2014-04-26 IHP
## 2014-04-27 automatically wrap mask, pop.args, det.args if not in list form
## 2014-09-03 linear mask tweaks
## 2014-11-23 groups
## 2014-11-24 new functions makeCH and processCH used by onesim()
## 2015-01-26 defaultextractfn updated
## 2015-01-26 streamlined outputtype (function getoutputtype called by both run.scenarios and fit.models)
## 2015-01-26 scen and repl argument for fit.models
## 2015-01-27 more robust handling of start values in processCH
## 2015-11-03 default extract fn reports unmarked and nonID sightings
## 2015-11-03 adapted for sim.resight etc.
## 2016-09-28 corrected error in working version that discarded supplied fit.args 'details'
## 2016-09-29 detectfn specified in det.args overrides scenario
## 2017-05-06 detectfn specified in scenario overrides det.args
## 2017-07-26 allow for S3 rbind.capthist
## 2017-09-14 nrepeat bug with method = 'none'
## 2017-12-01 tweak to fitarg$start: only one value per par
## 2020-01-28 change to new rse in defaultextractfn
## 2022-01-20 2.6.0 ncores default NULL; multithreaded secr.fit
## 2022-10-18 trapset components may be function; trap.args argument
## 2022-10-21 general tidy up
## 2022-12-28 allow fit.function = "ipsecr.fit"
## 2023-04-19 explicit fit = "multifit"
## 2023-04-29 maskset could be ignored in fitarg
## 2023-05-26 byscenario = TRUE fixed
## 2023-05-28 dynamic maskset for trapset function UNTESTED
## 2024-03-01 joinsessions argument
## 2024-05-01 is.function(trapset) messages
## 2024-09-27 pop.args model2D requires core
## 2025-06-07 add simOU.capthist option to args
## 2025-06-09 generalised detection parameters for OU

###############################################################################
wrapifneeded <- function (args, default) {
    if (any(names(args) %in% names(default)))
        list(args)  ## assume single naked list; wrap
    else
        args
}
###############################################################################
## complete a partial argument list (arg) with default args
## always return a list of arg lists long enough to match max(index)
fullargs <- function (args, default, index, multifit) {
    if (is.null(args)) {
        full.args <- list(default)
    }
    else {
        if (multifit) {
            # special case for multifit fit.arg 2023-04-14
            if (!is.list(args[[1]][[1]])) stop ("multifit expects nested list of fit.args")
            full.args <- mapply(fullargs, args, index = sapply(args,length), 
                MoreArgs = list(default = default, multifit = FALSE), SIMPLIFY = FALSE)
        }
        else {
            nind <- max(index)
            if (length(args) < nind) stop("too few components in args")
            tmpargs <- vector('list', nind)
            for (i in 1:nind) {
                if (is.character(args[[i]])) {
                    args[[i]] <- match.arg(args[[i]], default[[i]])
                }
                ## naked fn gives trouble here... 2014-09-03
                tmpargs[[i]] <- replace (default, names(args[[i]]), args[[i]])
                if (is.character(tmpargs[[i]]))
                    tmpargs[[i]] <- tmpargs[[i]][1]
            }
            full.args <- tmpargs
        }
    }
    full.args
}
###############################################################################

# 2023-02-07
designextractfn <- function(CH, ...) {
    if (!inherits(CH, 'capthist')) stop ("designextractfn expects capthist object")
    n <- nrow(CH)
    r <- sum(CH) - n
    dots <- list(...)
    esa <- sum(pdot(traps = traps(CH), ...)) * attr(dots$X, 'area')
    c(n = n, r = r, esa = esa, D = n/esa)
}
###############################################################################

defaultextractfn <- function(x, ...) {
    counts <- function(CH) {
        ## for single-session CH
        if (nrow(CH)==0) { ## 2015-01-24
            if (sighting(traps(CH)))
                c(n = 0, ndet = 0, nmov = 0, dpa = 0,
                  unmarked=0, nonID = 0, nzero = 0)
            else
                c(n=0, ndet=0, nmov=0, dpa = NA, rse = NA, rpsv = NA)
        }
        else {
            n <- nrow(CH)
            ndet <- sum(abs(CH)>0)
            r2 <- sum(abs(CH)) - n   ## 2020-01-28
            nmoves <- sum(unlist(sapply(moves(CH), function(y) y>0)))
            ## detectors per animal
            dpa <- if (length(dim(CH)) == 2)
                mean(apply(abs(CH), 1, function(y) length(unique(y[y>0]))))
            else
                mean(apply(apply(abs(CH), c(1,3), sum)>0, 1, sum))
            if (sighting(traps(CH))) {
                unmarked <- if (is.null(Tu <- Tu(CH))) NA else sum(Tu)
                nonID <- if (is.null(Tm <- Tm(CH))) NA else sum(Tm)
                nzero <- sum(apply(abs(CH),1,sum) == 0)
                c(n = n, ndet = ndet, nmov = nmoves, dpa = dpa,
                  unmarked=unmarked, nonID = nonID, nzero = nzero)
            }
            else {
                c(n=n, r=r2, nmov=nmoves, dpa = dpa, rse = 1 / sqrt(min(n,r2)), rpsv = RPSV(CH, CC = TRUE))
            }
        }
    }
    if (inherits(x, 'try-error')) {
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
    }
    else if (inherits(x, 'capthist')) {
        ## summarised raw data
        if (ms(x))
            unlist(lapply(x, counts))
        else {
            gp <- covariates(x)$group
            if (is.null(gp))
                counts(x)
            else
                unlist(lapply(split(x,gp,dropnullocc=TRUE), counts))
        }
    }
    else if (
        (inherits(x, 'secr') && !is.null(x$fit)) ||
        (inherits(x, 'ipsecr') && x$code == 1)
    ) {
        ## fitted model:
        ## default predictions of 'real' parameters
        out <- predict(x)
        if (!is.data.frame(out)) {
            warning ("summarising only first session, group or mixture class")
            out <- out[[1]]
        }
        attr(out, 'counts') <- counts(x$capthist)
        out
    }
    else {
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
    }
}

#####################
makeCH <- function (scenario, trapset, full.pop.args, full.det.args, 
    mask, multisession, joinsessions, detfunction) {
    ns <- nrow(scenario)
    # with( scenario, {
        CH <- vector(mode = 'list', ns)
        for (i in 1:ns) {
            #####################
            ## retrieve data
            grid   <- trapset[[scenario$trapsindex[i]]]
            poparg <- full.pop.args[[scenario$popindex[i]]]
            detarg <- full.det.args[[scenario$detindex[i]]]

            #####################
            ## POPULATION
            ## override D, core, buffer
            if (inherits(mask, 'linearmask'))               ## force to linear...
                poparg$model2D <- 'linear'
            if (poparg$model2D %in% c('IHP', 'linear')) {   ## linear
                ## 2017-10-03 to allow user to specify core directly,
                ## make this assignment conditional
                if (!inherits(poparg$core, 'mask')) {
                    poparg$core <- mask
                }

                ## for 'linear' case we may want a constant numeric value
                if (!is.character(poparg$D) && !is.function(poparg$D) && 
                        (length(poparg$D)<nrow(mask))) {
                    poparg$D <- scenario$D[i]
                }
                if (is.function(poparg$D) && packageVersion('secr') < '4.5.8') {
                    stop("density function requires secr version >= 4.5.8")
                }
                if (scenario$nrepeats[i]!=1)
                    stop("nrepeats > 1 not allowed for IHP, linear")
            }
            else {
                # if (!inherits(poparg$core, 'mask')) {    # conditional 2023-11-07
                    poparg$core <- attr(mask, "boundingbox")
                # }
                poparg$D <- scenario$D[i] * scenario$nrepeats[i]  ## optionally simulate inflated density
            }
            poparg$buffer <- 0

            #####################
            ## generate population
            pop <- do.call(sim.popn, poparg)
            
            #####################
            ## CAPTHIST
            ## form dp for sim.capthist or sim.resight
            ## form par for starting values in secr.fit()
            ## 'par' does not allow for varying link or any non-null model (b, T etc.)
            pnames <-  parnames(scenario$detectfn[i])
            dp <- c(as.list(scenario[i,pnames]), recapfactor = scenario$recapfactor[i])
            if ('detectpar' %in% names(detarg) && !is.symbol(detarg$detectpar)) {
                dp <- replace (dp, names(detarg$detectpar), detarg$detectpar)
            }
            
            #####################
            ## override det args as required
            detarg$traps      <- grid
            detarg$popn       <- pop
            detarg$detectfn   <- scenario$detectfn[i]
            detarg$noccasions <- scenario$noccasions[i]
            if ("detectpar" %in% names(detarg))
                detarg$detectpar  <- dp
            else
                detarg$detparmat  <- dp
            
            ####################
            ## function-specific kludges
            
            # force sighting to secr::sim.resight
            if (sighting(grid)) {
                detfunction <- "sim.resight"   
            }
            
            if (detfunction == "simCH") {
                CHfun <- ipsecr::simCH
            }
            else {
                CHfun <- get(detfunction, envir = sys.frame())
            }
            if (detfunction == "sim.resight") {
                if (!is.null(markocc(grid))) {
                    detarg$noccasions <- length(markocc(grid))
                    if (detarg$noccasions != scenario$noccasions[i])
                        warning("length of markocc attribute overrides noccasions in scenario")
                }
            }
            
            #####################
            ## simulate detection
            CHi <- do.call(CHfun, detarg)
            if (joinsessions && ms(CHi)) CHi <- join(CHi)   ## 2024-03-01
            
            #####################
            
            if (!is.na(scenario$nrepeats[i]))
                attr(CHi, "n.mash") <- rep(NA, scenario$nrepeats[i])
           
            ## 2022-11-24
            ## remember this realisation of D from function
            attr(CHi, 'D') <- attr(detarg$popn, 'D')
            ##
            
            ## 2023-02-06
            ## remember detection parameters
            attr(CHi, 'detectpar') <- detarg$detectpar
            attr(CHi, 'detparmat') <- detarg$detparmat
            ## 
            
            CH[[i]] <- CHi
        }
        if (ns > 1) {
            ## assume a 'group' column is present if ns>1
            names(CH) <- 1:ns  # group?
            if (is.function(multisession)) {
                # CH <- multisession(CH, group)   # drop group 2024-05-21
                CH <- multisession(CH)
            }
            else if (multisession) {
                CH <- MS.capthist(CH)
                if (!is.null(scenario$group)) session(CH) <- scenario$group
                CH
            }
            else {
                nc <- sapply(CH, nrow)
                CH$verify <- FALSE
                ####################################
                class(CH) <- c("capthist", "list")   
                CH <- do.call(rbind, CH)
                ####################################
                covariates(CH)$group <- rep(scenario$group, nc)
                CH
            }
        }
        else {
            CH[[1]]
        }
    # })
}
#####################
processCH <- function (scenario, CH, fitarg, extractfn, fit, fitfunction, byscenario, ...) {
    if (fit == 'design') {
        if (nrow(scenario)>1) warning("ignoring multiple groups for 'design' option")
        extractfn(CH, 
            X = fitarg[[1]]$mask, 
            detectfn = scenario$detectfn[1],
            detectpar = attr(CH, 'detectpar'),
            noccasions = scenario$noccasions[1], ...)
    }
    else if (fit == 'multifit') {
        fits <- lapply(fitarg, processCH, 
            scenario = scenario, 
            CH = CH, 
            extractfn = identity, 
            fit = TRUE, 
            fitfunction = fitfunction, 
            byscenario = byscenario)   # do not pass ...
        if (fitfunction == 'secr.fit') {
            fits <- secrlist(fits)
            names(fits) <- sapply(fitarg, '[[', 'model')
        }
        else {
            warning('multifit for non-secr fit returns list of fits rather than secrlist')
        }
        extractfn(fits, ...)
    }
    else if (is.logical(fit) && !fit) {
        extractfn(CH, ...)
    }
    else {
        ## form par for starting values in secr.fit()
        ## 'par' does not allow for varying link or any non-null model (b, T etc.)
        ## D, lambda0, g0, sigma are columns in 'scenario'
        par <- with(scenario, {
            ## 2.9.2 test for numeric to avoid bad start values for D (below)
            if (!is.null(attr(CH, 'D')) && is.numeric(attr(CH, 'D'))) D <- mean(attr(CH, 'D'))
            wt <- D/sum(D)
            if (detectfn[1] %in% 14:19) {
                list(D = sum(D) * nrepeats, lambda0 = sum(lambda0*wt), sigma = sum(sigma*wt))
            }
            else if (detectfn[1] %in% 20) {
                # ad hoc 2025-06-09
                list(D = sum(D) * nrepeats, epsilon = sum(epsilon*wt), sigma = sum(sigma*wt), tau = sum(tau*wt))
            }
            else {
                list(D = sum(D) * nrepeats, g0 = sum(g0*wt), sigma = sum(sigma*wt))
            }
        })
        ## prepare arguments for secr.fit() or ipsecr.fit()
        fitarg$capthist <- CH

        if (byscenario) fitarg$ncores <- 1L
        if (is.null(fitarg$model)) {
            fitarg$model <- defaultmodel(fitarg$CL, fitarg$detectfn)
        }

        if (!is.null(fitarg$start) && (fitarg$start[1] == 'true')) {
            ## check to see if simple 'true' values will work
            ## requires intercept-only model for all parameters
            model <- eval(fitarg$model)
            if (!is.list(model)) model <- list(model)
            vars <- unlist(lapply(lapply(model, terms), attr, 'term.labels'))
            if (fitfunction == "secr.fit") {
                if (fitarg$CL) par$D <- NULL
                if ((length(vars) != 0) && (fitarg$method == 'none')) {
                    ## not yet ready for interspersed beta coef
                    stop("method = 'none' requires full start vector for complex models")
                }
                if ('h2' %in% vars) par$pmix <- 0.5
            }
            fitarg$start <- lapply(par, '[', 1)   ## 2017-12-01 only first value
        }
        ##-------------------------------------------------------------------
        ## 2015-11-03 code for overdispersion adjustment of mark-resight data
        ## no simulations in first iteration; defer hessian
        if (fitfunction == "secr.fit") {
            chatnsim <- fitarg$details$nsim
            if (abs(chatnsim)>0) {
                fitarg$details <- as.list(replace(fitarg$details, 'nsim', 0))
                fitarg$details <- as.list(replace(fitarg$details, 'hessian', FALSE))
            }
        }
        ##-------------------------------------------------------------------
        if (fitfunction == "secr.fit") {
            fit <- try(do.call(secr.fit, fitarg))
        }
        else {
            fit <- try(do.call(ipsecr::ipsecr.fit, fitarg))
        }
        ##-------------------------------------------------------------------
        ## code for overdispersion adjustment of mark-resight data
        if (!ms(CH) && sighting(traps(CH)) && !inherits(fit, 'try-error')) {
            if ((abs(chatnsim) > 0) &  (logLik(fit)>-1e9)) {
                fitarg$details$nsim <- abs(chatnsim)
                fitarg$details$hessian <- TRUE
                fit$call <- NULL
                fitarg$start <- fit
                if (chatnsim<0)
                    fitarg$method <- "none"
                fit <- try(do.call(fitfunction, fitarg))
            }
        }
        ##-------------------------------------------------------------------
        extractfn(fit, ...)
    }
}
#####################

getoutputtype <- function (output) {
    for (i in 1:length(output)) {
        typical <- output[[i]][[1]]  ## ith scenario, first replicate
        if (length(typical) > 0) break
    }
    if (length(typical) == 0) stop ("no results found")

    ## outputtype: secrfit, predicted, coef, numeric
    outputtype <-
        if (inherits(typical, 'secr'))
            'secrfit'
    else if (inherits(typical, 'secrlist'))
        'multifit'
    else if (inherits(typical, 'ipsecr'))
        'ipsecrfit'
    else if (inherits(typical, 'summary.secr'))
            'secrsummary'
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
        else
            'user'
    outputtype
}
#####################

getoutputclass <- function (outputtype) {
    switch (outputtype,
        secrfit     = c("fittedmodels", 'secrdesign', 'list'),
        ipsecrfit   = c("fittedmodels", 'secrdesign', 'list'),
        predicted   = c("estimatetables", 'secrdesign', 'list'),
        derived     = c("estimatetables", 'secrdesign', 'list'),
        regionN     = c("estimatetables", 'secrdesign', 'list'),
        coef        = c("estimatetables", 'secrdesign', 'list'),
        user        = c("estimatetables", 'secrdesign', 'list'),
        secrsummary = c("summary", 'secrdesign', 'list'),
        capthist    = c("rawdata", 'secrdesign', 'list'),
        selectedstatistics = c("selectedstatistics", 'secrdesign', 'list'),
        "list"      ## otherwise as character 2023-04-14
    )
}

###############################################################################
run.scenarios <- function (
    nrepl,  
    scenarios, 
    trapset, 
    maskset, 
    xsigma = 4,
    nx = 32, 
    pop.args, 
    CH.function = c("sim.capthist", "simOU.capthist", "simCH"),
    det.args, 
    fit = FALSE, 
    fit.function = c("secr.fit", "ipsecr.fit"),
    fit.args, 
    chatnsim = 0, 
    extractfn = NULL, 
    multisession = FALSE, 
    joinsessions = FALSE,
    ncores = NULL, 
    byscenario = FALSE, 
    seed = 123,  
    trap.args = NULL,
    prefix = NULL,
    ...) {

    #--------------------------------------------------------------------------
    onesim <- function (r, scenario) {
        ## only one mask an fitarg allowed per scenario
        fitarg <- full.fit.args[[scenario$fitindex[1]]]
        if (is.function(trapset[[1]])) {
            # create each detector layout for this simulation
            # trapset <- mapply (do.call, trapset, trap.args, SIMPLIFY = FALSE)
            for (i in scenario$trapsindex) {
                # replace only within scope of this function May 2024
                trapset[[i]] <- do.call(trapset[[i]], trap.args[[i]])
            }
            ## allow dynamic mask
            if (is.null(maskset)) {
                warning ('dynamic mask under development')
                maskset <- lapply(trapset, make.mask, buffer = xsigma * scenario$sigma[1], 
                                  nx = nx, type = 'trapbuffer')
            }
        }
        msk <- findarg(fitarg, 'mask', 1, maskset[[scenario$maskindex[1]]])
        if (fit == "multifit") {
            for (i in 1:length(fitarg))
                fitarg[[i]]$mask <- findarg(fitarg[[i]], 'mask', 1, maskset[[scenario$maskindex[1]]])
        }
        else {
            fitarg$mask <- msk
        }
        CH <- makeCH(scenario, trapset, full.pop.args, full.det.args,
                     msk, multisession, joinsessions, CH.function)
        processCH(scenario, CH, fitarg, extractfn, fit, fit.function, byscenario, ...)
    }
    #--------------------------------------------------------------------------
    runscenario <- function(x) {
        out <- lapply(1:nrepl, onesim, scenario = x)
        scenarionumber <- x$scenario[1]
        message("Completed scenario ", scenarionumber)
        if (!is.null(prefix) && !byscenario) {
            outputlist <- list(out)
            names(outputlist) <- as.character(scenarionumber) 
            saveRDS(makeoutput(outputlist, x), 
                    file = paste0(prefix, scenarionumber, '.RDS'))
        }
        out
    }
    ##--------------------------------------------------------------------------
    makeoutput <- function (output, scen) {
        outputtype <- getoutputtype(output)
        if (outputtype == 'selectedstatistics')
            ## collapse replicates within a scenario into a matrix
            output <- lapply(output, do.call, what = rbind)
        message("Completed in ", round((proc.time() - ptm)[3]/60,3), " minutes")
        desc <- packageDescription("secrdesign")  ## for version number
        value <- list (
            call      = cl,
            version   = paste('secrdesign', desc$Version),
            starttime = starttime,
            proctime  = (proc.time() - ptm)[3],
            scenarios = scen,
            trapset   = trapset,
            trap.args = trap.args,
            maskset   = if (is.null(uts)) maskset else NULL,
            xsigma    = xsigma,
            nx        = nx,
            pop.args  = pop.args,
            CH.function = CH.function,
            det.args  = det.args,
            fit       = fit,
            fit.function = fit.function,
            fit.args  = fit.args,
            extractfn = extractfn,
            seed      = seed,
            chatnsim  = chatnsim,
            nrepl     = nrepl,
            output    = output,
            outputtype = outputtype
        )
        class(value) <- getoutputclass(outputtype)
        if (outputtype == 'regionN')
            attr(value, 'regionsize') <- sapply(output, function(x) attr(x[[1]], 'regionsize'))
        value
    }
    ##--------------------------------------------------------------------------
    
    ## mainline
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    
    # not forcing match.arg for CH.function allows user function
    CH.function <- CH.function[1]
    fit.function <- match.arg(fit.function)
    if (any(scenarios$detectfn == 20) && CH.function != "simOU.capthist") {
        CH.function <- "simOU.capthist"
        warning ("Using CH.function 'simOU.capthist' for Ornstein-Uhlenbeck detectfn 20")
    }
    if (CH.function == "simOU.capthist" && !any(scenarios$detectfn == 20)) {
        stop ("CH.function 'simOU.capthist' requires scenarios$detectfn 20")
    }
        starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    ncores <- secr::setNumThreads(ncores)   ## 2022-12-29
    if (byscenario && (ncores > nrow(scenarios)))
        stop ("when allocating by scenario, ncores should not exceed number of scenarios")
    if ((is.function(multisession) || multisession) && !anyDuplicated(scenarios$scenario)) {
        warning ("multisession ignored because no scenario duplicated")
    }

    ##--------------------------------------------
    ## default extractfn
    if (is.null(extractfn)) {
        if (fit == 'design')
            extractfn <- designextractfn
        else
            extractfn <- defaultextractfn
    }
    ##--------------------------------------------
    ## preprocess inputs
    if (inherits(trapset, 'traps') || is.function(trapset))  ## otherwise assume already list of traps
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

    if (any(sapply(trapset, is.function))) {
        if (!all(sapply(trapset, is.function)))
            stop ("all trapset must be a function if any is a function")
        if (missing(maskset)) 
            stop ("maskset should be provided if trapset is list of functions")
        if (length(trapset) != length(trap.args)) {
            stop ("trapset is list of functions, trap.args should be a list of the same length")
        }
        temptrapset <- list()
        for (i in 1:length(trapset)) {
            message ("Testing trapset function ", i, "...")
            temptrapset[[i]] <- do.call(trapset[[i]], trap.args[[i]])
        }
        dettype <- sapply(temptrapset, detector)[scenarios$trapsindex]
        sight <- sapply(temptrapset, function(x)
            if(ms(x)) sighting(x[[1]]) else  sighting(x))
        message ("Testing complete")
        message ("Detectors ", paste(
            sapply(temptrapset, function(x) if(ms(x)) nrow(x[[1]]) else nrow(x)),
            collapse = ', '))
    }
    else {
        dettype <- sapply(trapset, detector)[scenarios$trapsindex]
        sight <- sapply(trapset, function(x)
            if(ms(x)) sighting(x[[1]]) else  sighting(x))
    }
    
    if (!(all(sight) | all(!sight)))
        stop ("cannot mix sighting and nonsighting simulations")
    sight <- any(sight)

    OK <- !any((scenarios$nrepeats>1) & (dettype == "single"))
    OK <- if(is.na(OK)) TRUE else OK
    if (!OK)
        warning("single-catch traps violate independence assumption for nrepeats > 1")

    ##---------------------------------------------
    ## POPULATION ARGS
    ## allow user changes to default sim.popn arguments
    default.args <- as.list(args(sim.popn))[1:12]
    default.args$model2D  <- eval(default.args$model2D)[1]   ## 2014-09-03
    if (missing(pop.args)) pop.args <- NULL
    pop.args <- wrapifneeded(pop.args, default.args)
    full.pop.args <- fullargs (pop.args, default.args, scenarios$popindex, FALSE)
    
    ##---------------------------------------------
    ## CAPTHIST ARGS
    ## allow user changes to default sim.capthist or sim.resight arguments
    if (CH.function == "simCH") {
        if (!requireNamespace("ipsecr")) {
            stop ("requires package ipsecr; please install")
        }
        CHfun <- ipsecr::simCH
    }
    else {
        CHfun <- get(CH.function, envir = sys.frame())
    }
    default.args <- as.list(formals(CHfun))
    default.args[["..."]] <- NULL   # not relevant
    # if (sight)
    #     default.args <- as.list(args(sim.resight))[c(2,4:11)]  ## drop traps & ... argument
    # else
    #     default.args <- as.list(args(sim.capthist))[1:15]

    if (missing(det.args)) det.args <- NULL
    det.args <- wrapifneeded(det.args, default.args)
    full.det.args <- fullargs (det.args, default.args, scenarios$detindex, FALSE)

    ##---------------------------------------------
    ## FIT ARGS
    ## allow user changes to default fit.function arguments
    if (fit.function == 'secr.fit') {
        default.args <- as.list(formals(secr.fit))
        default.args[["..."]] <- NULL   # not relevant
        default.args$verify    <- FALSE   ## never check
        default.args$start     <- "true"  ## known values
        default.args$detectfn  <- 0       ## halfnormal
        default.args$biasLimit <- NA      ## never check
        default.args$details   <- list(nsim = 0)
        default.args$trace     <- FALSE
    }
    else if (fit.function == 'ipsecr.fit') {
        if (!requireNamespace("ipsecr")) stop ("requires package ipsecr; please install")
        default.args <- as.list(formals(ipsecr::ipsecr.fit))
        default.args[["..."]] <- NULL   # not relevant
        default.args$verify   <- FALSE    ## never check
        default.args$start    <- "true"   ## known values
        default.args$detectfn <- 0        ## halfnormal
        default.args$proxyfn  <- ipsecr::proxy.ms
        default.args$verbose  <- FALSE
    }
    else stop ("unrecognised fit function")
    if (missing(fit.args)) fit.args <- NULL
    fit.args <- wrapifneeded(fit.args, default.args)
    
    full.fit.args <- fullargs (fit.args, default.args, scenarios$fitindex, fit == "multifit")
    if (fit.function == "secr.fit") {
        for (i in 1:length(full.fit.args)) {
            if ('details' %in% names(full.fit.args[[i]]))
                full.fit.args[[i]]$details$nsim <- replace(full.fit.args$details,'nsim',chatnsim)
            ## stop("chatnsim not currently available for multifit models")
        }
    }
    
    ##--------------------------------------------
    ## MASKS
    ## construct masks as required
    if (missing(maskset)) {
        trapsigma <- scenarios[, c('trapsindex', 'sigma'), drop = FALSE]
        uts <- unique (trapsigma)
        code <- apply(uts, 1, paste, collapse='.')
        scenarios$maskindex <- match(apply(trapsigma,1,paste,collapse='.'), code)
        maskset <- vector('list', nrow(uts))
        for (k in 1:nrow(uts))
            maskset[[k]] <- make.mask(trapset[[uts$trapsindex[k]]],
                                      buffer = uts$sigma[k] * xsigma,
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
    ## check 'group'  (moved down to follow maskindex 2022-11-30)
    if ('group' %in% names(scenarios)) {
        scenarios$group <- factor(scenarios$group)
        if (any(tabulate(scenarios$group)>1)) {
            ## check no change of trapsindex etc. within group
            fields <- c('trapsindex','noccasions','nrepeats','fitindex','maskindex')
            fixed <- scenarios[,fields]
            scens <- split(fixed, scenarios$scenario)
            varyingfn <- function (x) apply(x, 2, function (y) length(unique(y))>1)
            varying <- sapply(scens, varyingfn)
            if (any(varying)) {
                cat ("Fields varying among groups within scenario - \n")
                print(varying)
                stop ("Fields ", paste(fields, collapse=', '), " must be constant across groups")
            }
        }
    }
    
    #--------------------------------------------
    ## override nrepeats and D in scenarios when IHP distribution
    for (i in 1:nrow(scenarios)) {
        pi <- scenarios$popindex[i]
        mi <- scenarios$maskindex[i]
        if ((full.pop.args[[pi]]$model2D %in% c('IHP', 'linear'))) {  ## linear 2014-09-03
            avD <- NA
            if (is.character(full.pop.args[[pi]]$D)) {          
                # avD <- mean (covariates(maskset[[mi]])[,full.pop.args[[pi]]$D])
                # bug 2024-05-19 does not have core at this point if core not in poparg
                # 2024-09-27 catch
                if (is.null(full.pop.args[[pi]]$core) || !inherits(full.pop.args[[pi]]$core, 'mask'))
                    stop ("pop.args: for model2D = 'IHP' with character 'D' specify a mask as argument 'core'")
                avD <- mean (covariates(full.pop.args[[pi]]$core)[,full.pop.args[[pi]]$D])
            }
            else if (is.numeric(full.pop.args[[pi]]$D)) {
                avD <- mean(full.pop.args[[pi]]$D)
            }
            scenarios[i, 'nrepeats'] <- 1   ## override
            if (!is.na(avD)) scenarios[i, 'D'] <- avD
        }
    }

    #--------------------------------------------
    ## run simulations
    tmpscenarios <- split(scenarios, scenarios$scenario)
    if (ncores > 1 && byscenario) {
        list(...)    ## ensures promises evaluated see parallel vignette 2015-02-02
        clustertype <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
        clust <- parallel::makeCluster(ncores, type = clustertype, methods = TRUE)
        if (clustertype == "PSOCK") {
            clusterEvalQ(clust, library(secr))
            clusterExport(clust, c(
                "runscenario", "onesim", "full.fit.args", "findarg",
                "maskset", "trapset", "trap.args", "full.det.args", 
                "multisession", "joinsessions", "CH.function", "makeCH", 
                "processCH", "extractfn", "fit", "fit.function", 
                "byscenario", ...
            ), environment())
        }
        parallel::clusterSetRNGStream(clust, seed)
        on.exit(parallel::stopCluster(clust))
        output <- parallel::parLapply(clust, tmpscenarios, runscenario)
    }
    else {
        set.seed (seed)
        output <- lapply(tmpscenarios, runscenario)
    }
    ##-------------------------------------------
    ## tidy output
    makeoutput (output, scenarios) 
    
}

########################################################################################

## version of run.scenarios that accepts existing data and
## expands scenarios for multiple model definitions

fit.models <- function (
    rawdata, 
    fit = FALSE, 
    fit.function = c("secr.fit", "ipsecr.fit"),
    fit.args, 
    chatnsim = 0, 
    extractfn = NULL,
    ncores = NULL, 
    byscenario = FALSE, 
    scen, 
    repl, 
    ...) {
    
    #--------------------------------------------------------------------------
    onesim <- function (CH, scenario) {
        fitarg <- full.fit.args[[scenario$fitindex[1]]]
        if (is.null(fitarg$mask)) {  
            fitarg$mask <- maskset[[scenario$maskindex[1]]]
        }
        processCH(scenario, CH, fitarg, extractfn, fit, fit.function, byscenario, ...)
    }
    #--------------------------------------------------------------------------
    runscenario <- function(x) {
        ## match by name, not number 2015-01-27
        scenID <- as.character(trunc(x$scenario[1]))
        out <- lapply(CHlist[[scenID]], onesim, scenario = x)
        message("Completed scenario ", x$scenario[1])
        flush.console()
        out
    }
    ##--------------------------------------------------------------------------
    ## mainline
    fit.function <- match.arg(fit.function)
    if (!inherits(rawdata, "rawdata"))
        stop ("requires rawdata output from run.scenarios()")
    ## optionally select which scenarios to fit
    if (missing(scen)) {
        scen <- unique(rawdata$scenarios$scenario)
    }

    else {
        scen <- unique(scen)
        if (!all(scen %in% unique(rawdata$scenarios$scenario)))
            stop ("invalid scen argument")
        if (length(scen)<1)
            stop ("invalid scen argument")
    }
    ## optionally select which replicates to fit
    if (missing(repl)) {
        nrepl <- rawdata$nrepl
        repl <- 1:nrepl
    }
    else {
        repl <- unique(repl)
        if (!all(repl %in% 1:rawdata$nrepl))
            stop ("invalid repl argument")
        nrepl <- length(repl)
        if (nrepl<1)
            stop ("invalid repl argument")
    }
    CHlist <- lapply(rawdata$output[scen], '[', repl)
    scenarios <- rawdata$scenarios[rawdata$scenarios$scenario %in% scen,]

    trapset <- rawdata$trapset
    maskset <- rawdata$maskset
    
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    # if (is.null(ncores)) {
    #     ncores <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
    # }
    ncores <- secr::setNumThreads(ncores)   ## 2022-12-29
    if (byscenario & (ncores > nrow(scenarios))) {
        stop ("ncores exceeds number of scenarios")
    }

    ##--------------------------------------------
    ## default extractfn
    if (is.null(extractfn)) {
        extractfn <- defaultextractfn
    }

    ##---------------------------------------------
    ## allow user changes to default arguments
    if (fit.function == 'secr.fit') {
        default.args <- as.list(args(secr.fit))[1:21]
        default.args$biasLimit <- NA       ## never check
        default.args$verify    <- FALSE    ## never check
        default.args$start     <- "true"   ## known values
        default.args$detectfn  <- 0        ## halfnormal
        default.args$details   <- list(nsim = 0)
        default.args$trace     <- FALSE
    }
    else if (fit.function == 'ipsecr.fit') {
        if (!requireNamespace("ipsecr")) stop ("requires package ipsecr; please install")
        default.args <- as.list(formals(ipsecr::ipsecr.fit))[1:16]
        default.args$proxyfn <- ipsecr::proxy.ms
        default.args$verify  <- FALSE   ## never check
        default.args$start   <- "true"  ## known values
        default.args$verbose <- FALSE
    }
    else stop ("unrecognised fit function")

    if (missing(fit.args)) fit.args <- NULL
    fit.args <- wrapifneeded(fit.args, default.args)
    nfit <- length(fit.args)
    if (nfit > 1) {
        ## expand scenarios by the required number of different model fits
        ##      scenarios <- scenarios[rep(scenarios$scenario, each = nfit),]
        scenarios <- scenarios[rep(1:nrow(scenarios), each = nfit),]
        scenarios$fitindex <- rep(1:nfit, length.out = nrow(scenarios))
        ## assign new unique scenario number by adding decimal fraction
        scenarios$scenario <- scenarios$scenario + scenarios$fitindex /
            10 ^ trunc(log10(nfit)+1)
        scenarios <- scenarios[order(scenarios$scenario),]
        rownames(scenarios) <- 1:nrow(scenarios)

    }
    full.fit.args <- fullargs (fit.args, default.args, scenarios$fitindex, fit == "multifit")

    for (i in 1: length(full.fit.args))
        full.fit.args[[i]]$details <- as.list(replace(full.fit.args[[i]]$details,'nsim',chatnsim))

    ## construct masks as required
    if (is.null(maskset)) {
        trapsigma <- scenarios[, c('trapsindex', 'sigma'), drop = FALSE]
        uts <- unique (trapsigma)
        code <- apply(uts, 1, paste, collapse='.')
        scenarios$maskindex <- match(apply(trapsigma,1,paste,collapse='.'), code)
        maskset <- vector('list', nrow(uts))
        for (k in 1:nrow(uts))
            maskset[[k]] <- make.mask(trapset[[uts$trapsindex[k]]], buffer = uts$sigma[k] *
                                      rawdata$xsigma, type = 'trapbuffer', nx = rawdata$nx)
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
    ## run simulations
    tmpscenarios <- split(scenarios, scenarios$scenario)

    if (ncores > 1 && byscenario) {
        list(...)    ## ensures promises evaluated see parallel vignette 2015-02-02
        clustertype <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
        clust <- parallel::makeCluster(ncores, type = clustertype, methods = TRUE)
        if (clustertype == "PSOCK") {
            clusterEvalQ(clust, library(secr))
            clusterExport(clust, c(
                "runscenario", "onesim", "full.fit.args", "findarg",
                "maskset", "trapset", "trap.args", "full.det.args", 
                "multisession", "joinsessions", "CH.function", "makeCH", 
                "processCH", "extractfn", "fit", "fit.function", 
                "byscenario", ...
            ), environment())
        }
        on.exit(parallel::stopCluster(clust))
        output <- parallel::parLapply(clust, tmpscenarios, runscenario)
    }
    else {
        output <- lapply(tmpscenarios, runscenario)
    }
    
    ##-------------------------------------------
    ## tidy output
    
    outputtype <- getoutputtype(output)
    if (outputtype == 'selectedstatistics')
        ## collapse replicates within a scenario into a matrix
        output <- lapply(output, do.call, what = rbind)
    message("Completed in ", round((proc.time() - ptm)[3]/60,3), " minutes")
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
                   nrepl = nrepl,     ## rawdata$nrepl, 2015-01-26
                   output = output,
                   outputtype = outputtype
                   )
    class(value) <- getoutputclass (outputtype)
    if (outputtype == 'regionN')
        attr(value, 'regionsize') <- sapply(output, function(x) attr(x[[1]], 'regionsize'))

    value
}
