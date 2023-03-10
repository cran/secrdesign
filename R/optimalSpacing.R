##############################################################################
## package 'secrdesign'
## optimalSpacing.R
## 2017-07-15 moved from Lambda.R
## 2017-08-09 finished (?) fiddling with plot etc.
## 2018-11-28 distribution; openCR
## 2019-02-15 try-error catch bad uniroot in interpRSE
## 2020-01-20 2.6.0 removed openCR
##############################################################################

interpRSE <- function (values) {
    ## find minimum at n = r
    nminr <- function(RR) {
        n <- approx(values$R, values$n, RR)$y
        r <- approx(values$R, values$r, RR)$y
        n-r
    }
    rangeR <- range(values$R)
    if (diff(rangeR)>0) {
        ur <- try(uniroot(nminr, interval = rangeR))
        if (inherits(ur, 'try-error'))
            list(minimum = NA, objective = NA)
        else
            list(minimum = ur$root, objective = approx(values$R, values$RSE, ur$root)$y)
    }
    else         list(minimum = NA, objective = NA)
    
}

oneRSE <- function (R, k, traps, xsigma, detectpar, noccasions, nrepeats, detectfn,
                                            maskargs = list(), CF) {
    # hold traps constant, vary sigma, D, mask
    tmpsigma <- spacing(traps) / R
    tmpD <- (k / tmpsigma)^2
    detectpar$sigma <- tmpsigma
    maskargs$buffer <- xsigma * tmpsigma
    maskargs$spacing <- maskargs$spacing / R
    mask <- do.call(make.mask, maskargs)
    if (nrow(mask) < 600) warning ("mask has only ", nrow(mask), " points")
    minnrRSE(D = tmpD * nrepeats, traps, mask, detectpar, noccasions, detectfn, CF)
}
##############################################################################

getrotRSE <- function (R, k, traps, xsigma, detectpar, noccasions, nrepeats,
                       detectfn, maskargs = list(), CF, 
                       distribution = "poisson") {

    tmpsigma <- spacing(traps) / R
    tmpD <- (k / tmpsigma)^2
    maskargs$buffer <- tmpsigma * xsigma
    maskargs$spacing <- maskargs$spacing / R
    mask <- do.call(make.mask, maskargs)
    detectpar$sigma <- tmpsigma

    nrm <- Enrm(tmpD, traps = traps, mask = mask, detectpar = detectpar,
                noccasions = noccasions, detectfn = detectfn) * nrepeats

    rotRSE <- minnrRSE(tmpD * nrepeats, traps = traps, mask = mask, detectpar = detectpar,
                        noccasions = noccasions, detectfn = detectfn, CF = 1.0)
    
    if (distribution == "binomial") {
        # 2019-01-08 blocked as doesn't seem to do anything
        # Pxy <- pdot(mask, traps, detectfn, detectpar, noccasions)
        # esa <- sum(Pxy) * attr(mask, 'area')
        # 2019-01-08 inserted * nrepeats 
        rotRSE <- sqrt(rotRSE^2 - 1 / (tmpD * nrepeats * maskarea(mask)))
    }
    c(R, nrm, rotRSE*CF)
}
##############################################################################

simRSEfn <- function (R, k, traps, xsigma, detectpar, noccasions, nrepeats, detectfn,
                      nx, nrepl, fitfunction, allargs) {
    tmpsigma <- spacing(traps) / R  # vector same length as R
    tmpD <- (k / tmpsigma)^2        # vector same length as R
    
    ## scenarios
    scen1 <- make.scenarios(trapsindex = 1, noccasions = noccasions, nrepeats = nrepeats, D = tmpD,
                            sigma = 1, lambda0 = detectpar$lambda0, detectfn = detectfn)
    scen1$sigma <- tmpsigma  ## specify here to avoid crossing with D
    
    ## general arguments
    defaultargs <- list(nrepl = nrepl, trapset = traps, scenarios = scen1, fit = TRUE,
                        fit.function = fitfunction, byscenario = FALSE, ncores = NULL, 
                        xsigma = xsigma, nx = nx)
    dotsargs <- allargs[names(allargs) %in% c("seed", "ncores")]
    runargs <- replacedefaults(defaultargs, dotsargs)

    ## sim.popn arguments: only Ndist
    defaultpopargs <- list(Ndist = "poisson")
    popargs <- allargs[names(allargs) %in% c("Ndist")]
    runargs$pop.args <- replacedefaults(defaultpopargs, popargs)
    
    ## arguments for model fitting function
    defaultfitargs <- list(start = "true", 
                           detectfn = detectfn, 
                           details = list())
    fitargs <- allargs[!(names(allargs) %in% c("seed", "ncores", "Ndist"))]
    fitargs <- fitargs[names(fitargs) %in% names(formals(fitfunction))]   ## 2019-02-15
    runargs$fit.args <- replacedefaults(defaultfitargs, fitargs)    
    ## run
    sims1 <- do.call(run.scenarios, runargs)
    
    allRSE <- select.stats(sims1, parameter = "D", c("estimate","RSE"))
    eachRSE <- cbind(R = rep(R, each = nrepl), do.call(rbind, allRSE$output))
    tmp <- summary(select.stats(sims1, parameter = "D", c("RSE","RB","ERR")),
                   fields = c("n","mean","se","rms"))$OUTPUT
    tmp2 <- lapply(tmp, unlist)
    simout <- as.data.frame(t(sapply(tmp2, '[', c(1,4,7,5,8,12))))
    names(simout) <- c("n","RSE.mean","RSE.se","RB.mean", "RB.se","rRMSE")
    simout$rRMSE <- simout$rRMSE / scen1$D
    list(eachRSE = eachRSE, summary = cbind(data.frame(R = R), simout))
}
##############################################################################

optimalSpacing <- function (
    D, 
    traps, 
    detectpar, 
    noccasions,
    nrepeats = 1,
    detectfn = c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX'),
    fittedmodel = NULL,
    xsigma = 4,
    R = seq(0.2, 4, 0.2),
    CF = 1.0,
    distribution = c("poisson", "binomial"),
    fit.function = c("none", "secr.fit"),
    simulationR = seq(0.4, 4, 0.4),
    nrepl = 10,
    plt = FALSE,
    ...) {

    ## if a fitted model is provided then all preceding arguments are overridden
    ## (D, traps, detectpar, noccasions, detectfn)

    ## This function seeks the detector spacing that minimises RSE for a given geometry.
    ## It uses the trick of
    ## (i)  expressing spacing as a multiple of sigma, and
    ## (ii) varying (relative) spacing by scaling sigma and other spatial parameters
    ##      while keeping the detector array fixed.
    ## For relative spacing = R, the scaling factors are
    ##    sigma = sigma / R
    ##    D = (k / sigma)^2
    ##    buffer = xsigma * sigma / R
    ##    maskspacing = maskspacing / R

    ## note that if trap spacing is doubled,
    ## - the number of animals per trap square is quadrupled
    ## - sigma / spacing is halved

    ## suppress simulations with nrepl = 0
    ## suppress optimisation with CF = NA
    ## control values with R
    
    fit.function <- match.arg(fit.function)
    distribution <- match.arg(distribution)
    if (!is.null(fittedmodel)) {
        if (ms(fittedmodel$capthist))
            stop ("optimalSpacing requires single-session data")
        pred <- predict(fittedmodel)
        detectfn <- fittedmodel$detectfn
        detectpar = list(g0 = pred['g0','estimate'],
                         lambda0 = pred['lambda0','estimate'],
                         sigma = pred['sigma','estimate'],
                         z = pred['z','estimate'],
                         w = pred['w','estimate'])

        traps <- attr(fittedmodel$capthist, 'traps')
        if (any(detector(traps)=='single'))
            warning ("results are unreliable for single-catch traps when lambda0 ",
                     "inferred from multi-catch model")
        noccasions <- ncol(fittedmodel$capthist)

        if (fittedmodel$CL) {
            D <- derived(fittedmodel)['D','estimate']
        }
        else {
            D <- pred['D','estimate']
        }
    }
    else {
        if (is.character(detectfn)) detectfn <- match.arg(detectfn)
        detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:19))

    }
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar

    args <- list(...)
    k <- detectpar$sigma * D^0.5

    # prepare make.mask arguments
    defaultmaskargs <- list (traps = traps, buffer = xsigma * spacing(traps),
                             nx = 64, type = "trapbuffer")
    dotsargs <- args[names(args) %in% c('nx', 'type', 'poly','poly.habitat')]
    maskargs <- replacedefaults(defaultmaskargs, dotsargs)
    tmpmask <- do.call(make.mask, maskargs) # baseline mask constructed at R = 1
    maskargs$spacing <- spacing(tmpmask)

    if (any(detector(traps) == 'single')) {
        warning ("treating single-catch traps as multi-catch", call. = FALSE)
        detector(traps) <- 'multi'
    }

    #################
    values <- sapply(R, getrotRSE, k, traps, xsigma, detectpar,
                     noccasions, nrepeats, detectfn, maskargs, CF, distribution)
    values <- as.data.frame(t(values))
    names(values) <- c("R", "n", "r", "m", "RSE")
    #################

    if (!is.na(CF)) {
        if (!is.null(values))
            opt <- try(interpRSE(values))
        else
            ## this path is never taken in current version 2017-09-27
            opt <- optimize(oneRSE,  interval = range(R), k = k, traps = traps,
                            xsigma = xsigma, detectpar = detectpar, noccasions = noccasions,
                            nrepeats = nrepeats, detectfn = detectfn, maskargs = maskargs,
                            CF = CF)
    }
    #################
    if (fit.function %in% c("secr.fit")) {
        args$details <- list(distribution = distribution)   ## blocks other supplied details
        simRSE <- simRSEfn (simulationR, k, traps, xsigma, detectpar, noccasions, nrepeats,
                            detectfn, nx = maskargs$nx, nrepl, fit.function, args)
    }
    else simRSE <- NULL
    #################

    rotRSE <- list(values = values)
    if (!is.na(CF)) {
        rotRSE$optimum.spacing <- opt$minimum * detectpar$sigma # spacing(traps)
        rotRSE$optimum.R <- opt$minimum
        rotRSE$minimum.RSE <- opt$objective
    }
    else {
        rotRSE$optimum.spacing <- rotRSE$optimum.R <- rotRSE$minimum.RSE <- NA
    }

    out <- list(rotRSE = rotRSE,
                simRSE = simRSE)
    attr(out, "noccasions") <- noccasions
    attr(out, "traps") <- traps
    attr(out, "nrepeats") <- nrepeats
    attr(out, "detectpar") <- detectpar
    attr(out, "detectfn") <- detectfn

    class(out) <- c("optimalSpacing", "list")

    if (plt) {
        ## fixed 2017-08-24
        args$x <- out
        if (is.null(args$add)) args$add <- FALSE
        if (is.null(args$plottype)) args$plottype <- "RSE"
        do.call(plot, args)
        invisible(out)
    }
    else {
        out
    }
}
##############################################################################

plot.optimalSpacing <- function (x, add = FALSE, plottype = c("RSE", "nrm"), ...) {
    ## need to define missing cases
    args <- list(...)
    plottype <- match.arg(plottype)
    if (plottype == "nrm") {
        y <- x$rotRSE$values$n + x$rotRSE$values$r
    }
    else {
        y <- x$rotRSE$values$RSE
        if(all(is.na(y))) {
            warning ("RSE all NA")
        }
    }

    R <- x$rotRSE$values$R
    if (!add) {
        maxy <- 0.5
        if (!all(is.na(y))) {
            maxy <- max(y, na.rm = TRUE)*1.3
        }
        maxx <- max(R, na.rm = TRUE)
        minx <- min(R, na.rm = TRUE)
        if (minx < 0.2 * (maxx-minx)) minx <- 0
        defaultargs <- list(x = 0, y = 0, type = "n", las = 1,
                            xlab = expression(paste("Spacing -  ", sigma, "  units")),
                            ylab = expression(paste("RSE ", hat(italic(D)))),
                            ylim = c(0, maxy),
                            xlim = c(minx, maxx))
        if (plottype == 'nrm') defaultargs$ylab <- "Number"
        dotsargs <- args[names(args) %in% c("xlab", "ylab", "xlim", "ylim", "las",
                                            "xaxs", "yaxs")]
        plotargs <- replacedefaults(defaultargs, dotsargs)
        do.call(plot, plotargs)
    }

    if (plottype == "RSE") {
        defaultargs <- list(col = "black", lwd = 1, cex = 1, pch = 16)
        dotsargs <- args[names(args) %in% c("col", "lwd", "lty", "cex", "pch")]
        plotargs <- replacedefaults(defaultargs, dotsargs)

        plotargs$x <- R
        plotargs$y <- y
        do.call(lines, plotargs)

        # suspend this 2018-11-30
        # if (!is.null(x$rotRSE$optimum.R)) {
        #     plotargs$x <- x$rotRSE$optimum.R
        #     plotargs$y <- x$rotRSE$minimum.RSE
        #     do.call(points, plotargs)
        # }
        if (!is.null(x$simRSE)) {
            plotargs$pch <- 1
            plotargs$x <- x$simRSE$summary$R
            plotargs$y <- x$simRSE$summary$RSE.mean
            ## 2017-09-18
            segments(plotargs$x, plotargs$y - 2 * x$simRSE$summary$RSE.se,
                     plotargs$x, plotargs$y + 2 * x$simRSE$summary$RSE.se)
            do.call(points, plotargs)
        }
    }
    if (plottype == "nrm") {
        defaultargs <- list(col = "blue", lwd = 1, cex = 1, pch = 16)
        dotsargs <- args[names(args) %in% c("col", "lwd", "lty", "cex", "pch")]
        plotargs <- replacedefaults(defaultargs, dotsargs)

        plotargs$x <- R
        plotargs$y <- x$rotRSE$values$n
        do.call(lines, plotargs)

        plotargs$col <- "red"
        plotargs$x <- R
        plotargs$y <- x$rotRSE$values$r
        do.call(lines, plotargs)

        plotargs$lty <- 2
        plotargs$x <- R
        plotargs$y <- x$rotRSE$values$m
        do.call(lines, plotargs)

        legend (x = "top", legend = c("n","r","m"), lty=c(1,1,2), col=c("blue","red","red"), horiz = TRUE,
                cex = par()$cex * 1.2)
    }

}

print.optimalSpacing <- function (x, ...) {
    attr(x,"class") <- NULL
    print(x)
}
