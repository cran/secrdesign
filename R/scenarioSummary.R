# Deterministic summary of scenarios

costing <- function(traps, nr, noccasions, unitcost = list(), nrepeats = 1, routelength = NULL, 
                    setupoccasion = TRUE) {

    defaultunitcost <- list(perkm = 0, perarray = 0, perdetector = 0, pervisit = 0, perdetection = 0)
    unitcost <- replacedefaults (defaultunitcost, unitcost)
    
    trapcost <- covariates(traps)$costpervisit
    if (is.null(trapcost)) {
        trapcost <- rep(1, nrow(traps))
    }
   
    if (is.null(routelength)) {
        spc <- spacing(traps)
        if (is.null(spc)) spc <- 0   ## bug fix 2019-01-09
        routelength <- spc * (nrow(traps)-1) / 1000
    }
    
    nocc <- noccasions + setupoccasion

    costs <- c(travel = unitcost$perkm * routelength * nocc * nrepeats,
               arrays = unitcost$perarray * nrepeats,
               detectors = unitcost$perdetector * nrow(traps) * nrepeats,
               visits = sum(unitcost$pervisit * trapcost) * nocc * nrepeats,
               detections = unitcost$perdetection * sum(nr[1:2]))
    c(costs, totalcost = sum(costs))
}

trapsperHR <- function(traps, r) {
    mask <- make.mask(traps, buffer = r)
    d <- edist(mask, traps)
    nt <- apply(d, 1, function(x) sum(x<=r))
    nt <- nt[nt>0]
    c(median = median(nt), max = max(nt))
}
# 
# trapsperHR (gr, circular.r(p = 0.95, detectfn = 'HHN', sigma = 20))
# trapsperHR (gr, circular.r(p = 0.95, detectfn = 'HEX', sigma = 20))

scenarioSummary <- function (scenarios, trapset, maskset, xsigma = 4, nx = 64,   
                             CF = 1.0, costing = FALSE, ..., ncores = 1) {
    ## mainline
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    extrafields <- FALSE
    
    if (!all(as.character(scenarios$detectfn) %in% 
             c(as.character(c(0,1,2,14:18)), 'HN','HR','EX','HHN', 'HEX', 'HHR', 'HCG','HAN')))
        stop ("scenarioSummary requires hazard detection function HHN, HEX etc. (HN, EX approximated)")
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
    
    ##--------------------------------------------
    ## construct masks as required
    
    if (missing(maskset)) {
        trapsigma <- scenarios[, c('trapsindex', 'sigma'), drop = FALSE]
        uts <- unique (trapsigma)
        code <- apply(uts, 1, paste, collapse='.')
        scenarios$maskindex <- match(apply(trapsigma,1,paste,collapse='.'), code)
        maskset <- vector('list', nrow(uts))
        for (k in 1:nrow(uts)) {
            maskset[[k]] <- make.mask(trapset[[uts$trapsindex[k]]], 
                                      buffer = uts$sigma[k] * xsigma,
                                      type = 'trapbuffer', nx = nx)
        }
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
        
    }
    if (max(scenarios$maskindex) > length(maskset))
        stop ("maskindex does not match maskset")
    
    ##---------------------------------------------------------------------------
    ##
    detperHR <- function (traps, mask, detectfn, detectpar, noccasions, p=0.95) {
        pd <- pdot(mask, traps, detectfn, detectpar, noccasions)
        rad <- circular.r(p = p, detectfn = detectfn, detectpar = detectpar)
        ntxy <- function(xy, traps, rad) {
            sum (edist(traps, matrix(xy, nrow=1)) < rad)
        }
        nt  <- apply(mask, 1, ntxy, traps = traps, rad = rad)
        sumpd <- sum(pd)
        Ent <- sum(pd * nt) / sumpd
        Ent2 <- sum(pd * nt^2) / sumpd
        varnt <- Ent2 - Ent^2
        
        A <- pi * rad ^2
        ntratio <- sum(pd * nt/(A/spacing(traps)^2)) /sumpd
            
        Epd <- sum(pd^2) / sumpd
        Epd2 <- sum(pd^3) /sumpd
        varpd <- Epd2 - Epd^2
        
        # orient <- function (i,j) {
        #     diameter2 <- (detectpar$sigma * 2)^2
        #     dx <- traps$x[i]-traps$x[j]
        #     dy <- traps$y[i]-traps$y[j]
        #     out <- atan2 (dy, dx)
        #     out[(dx^2 + dy^2) > diameter2] <- NA
        #     out
        # } 
        # grad <- outer(1:nrow(traps), 1:nrow(traps), FUN = orient) 
        # SDorient <- sd(grad, na.rm = TRUE)
        # c(meannt = Ent, CVnt = sqrt(varnt)/Ent, SDorient = SDorient)
        c(meannt = Ent, CVnt = sqrt(varnt)/Ent, ntratio = ntratio, CVPxy = sqrt(varpd)/Epd)
    }
    
    saturation <- function (D, detectpar, detectfn, tr, mask) {
        if (!(detector(tr)[1] %in% c('multi','proximity', 'capped', 'count'))) {
            warning ("saturation only available for 'multi','proximity', 'capped' or 'count' detectors")
            return(NA)
        }
        if (is.character(detectfn))
            detectfn <- match.arg(detectfn)
        detectfn <- secr:::valid.detectfn(detectfn, 14:18)
        cellarea <- attr(mask, 'area') 
        dk <- edist(tr, mask)                                             # K x M
        if (detectfn == 14)
            lambda <- exp(-dk^2/2/detectpar$sigma^2)  # K x M
        else if (detectfn == 15)
            lambda <- 1 - exp(-(dk/detectpar$sigma)^-detectpar$z)
        else if (detectfn == 16)
            lambda <- exp(-dk^2/2/detectpar$sigma^2)  
        else if (detectfn == 17)
            lambda = exp(- (dk-detectpar$w)^2 / 2 / detectpar$sigma^2)
        else if (detectfn == 18)
            lambda = pgamma(dk, shape = detectpar$z, scale = detectpar$sigma/detectpar$z, 
                            lower.tail = FALSE, log.p = FALSE)
        lambda <- detectpar$lambda0 * lambda
        if (detector(tr)[1] == "multi") {
            Hi <- apply(lambda, 2, sum)                                 # M       overall hazard of detn animal at x (summed over traps)
            hmult <- (1 - exp(-Hi)) / Hi
            pkx <- sweep(lambda, MARGIN = 2, STATS = hmult, FUN = "*")  # K x M   Pr animal at x caught trap k
            Hk <- apply(-log(1-pkx), 1, sum) * D * cellarea             # K     integrated hazard for trap k. Assumes homogeneous D
            mean(1-exp(-Hk))  
        }
        else {
            H <- apply(lambda, 1, sum) * D * cellarea                         # K
            mean(1-exp(-H))
        }
        
    }
    
    ##---------------------------------------------------------------------------
    ##
    onescenario <- function (scenario) {

        traps <- trapset[[scenario$trapsindex]]
        mask <- maskset[[scenario$maskindex]]
        
        detectfn <- secr:::valid.detectfn (scenario$detectfn, c(0,1,2,14:18))
        detectpar = list(g0 = scenario$g0, lambda0 = scenario$lambda0, sigma = scenario$sigma, z = scenario$z)
        dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 1 to 15, 2 to 16
        detectfn <- dfc$detectfn
        detectpar <- dfc$detectpar
        
        Pxy <- pdot(mask, traps, detectfn, detectpar, scenario$noccasions)
        esa <- sum(Pxy) * attr(mask, 'area')
        
        sat <- saturation (scenario$D, detectpar, detectfn, traps, mask)
        if (detector(traps)[1] %in% c("multi", "proximity", "count"))
            nrm <- Enrm(scenario$D, traps, mask, detectpar, 
                        scenario$noccasions, detectfn = detectfn)
        else
            nrm <- c(En=NA,Er=NA,Em=NA)
        
        ## 2018-06-11
        nrm <- nrm * scenario$nrepeats
        esa <- esa * scenario$nrepeats
        
        rotRSE <- 1 / sqrt(min(nrm[1:2]))
        rotRSEB <- (rotRSE^2 - 1 / (scenario$D * scenario$nrepeats * maskarea(mask)))^0.5 * scenario$CF
        rotRSE <- rotRSE * scenario$CF

        out <- c(scenario[1:8], round(nrm,3))

        out <- c(out, 
                 esa = esa,
                 CF = scenario$CF, 
                 rotRSE = round(rotRSE,4),
                 rotRSEB = round(rotRSEB,4),
                 arrayN = ntraps[scenario$trapsindex],
                 arrayspace = round(spaces[scenario$trapsindex]/scenario$sigma,4),
                 arrayspan = round(spans[scenario$trapsindex]/scenario$sigma,4),
                 saturation = sat)

        if (costing) {
            dots <- list(...)
            if (is.null(dots$routelength)) 
                routelength <- NULL
            else { # extract routelength for current scenario
                if (length(dots$routelength)==1)
                    routelength <- dots$routelength
                else 
                    routelength <- dots$routelength[scenario$trapsindex]
                dots$routelength <- NULL  # drop from dots
            }
              
            arg <- list (traps = traps, 
                         nr = nrm, 
                         noccasions = scenario$noccasions,
                         nrepeats = scenario$nrepeats,
                         routelength = routelength)
            arg <- c(arg, dots)
            out <- c(out, do.call('costing', arg))
        }

        radius <- circular.r(p = 0.95, detectfn = detectfn, sigma = detectpar$sigma)
        out$detperHR <- trapsperHR (traps, radius)[1]  ## median detectors per 95% HR
        
        if (extrafields) {
            out <- c(out,
                     detperHR(traps, mask, detectfn = scenario$detectfn, detectpar, 
                              scenario$noccasions),
                     sinuosity = (spans/routelength)[scenario$trapsindex])
        }
        
        unlist(out)
    }
    ##---------------------------------------------------------------------------

    ## tweaked 2019-01-09 for single-trap arrays
    getspan <- function (traps) suppressWarnings(pmax(0, max(dist(traps))))
    spans <- unname(sapply(trapset, getspan))
    ntraps <- unname(sapply(trapset, nrow))
    spaces <- unname(sapply(trapset, spacing))
    spaces <- sapply(spaces, function(x) if (is.null(x)) NA else x)  ## patch 2019-01-09
    
    if (is.null(scenarios$CF)) {
        if (length(CF) > length(trapset))
            stop ("length of CF exceeds number of detector layouts")
        if (length(CF) < length(trapset))
            CF <- rep(CF, length(trapset))[1:length(trapset)]
        scenarios$CF <- CF[scenarios$trapsindex]
    }
    
    #--------------------------------------------
    ## run summaries
    tmpscenarios <- split(scenarios, scenarios$scenario)
    if (ncores > 1) {
        ## not needed: list(...) ## ensures promises evaluated see parallel vignette
        clust <- makeCluster(ncores, methods = TRUE)
        on.exit(stopCluster(clust))
        output <- parLapply(clust, tmpscenarios, onescenario)
    }
    else {
        output <- lapply(tmpscenarios, onescenario)
    }
    output <- do.call(rbind, output)
    output <- as.data.frame(output)
    output
}