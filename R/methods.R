###############################################################################
## package 'secrdesign'
## methods.R
## 2014-02-07,08
## 2014-04-09, 2014--04-11, 2014-04-29
###############################################################################

select.stats <- function (object, parameter = 'D', statistics) {
    if (!inherits(object, 'estimatetables'))
        stop ("select.stats requires input of class estimatetables")
    if (is.na(object$outputtype))
        stop ("cannot select.stats output of unknown type")
    if (object$outputtype == 'secrfit')
        stop ("cannot select.stats secr fit - use predict() first")

    estname <- ""
    SEname <- ""
    if (object$outputtype %in% c('predicted', 'derived', 'regionN')) {
        estname <- 'estimate'
        SEname = 'SE.estimate'
    }
    else if (object$outputtype == 'coef') {
        estname <- 'beta'
        SEname <- 'SE.beta'
    }

    typical <- object$output[[1]][[1]]  ## first scenario, first replicate
    stat0 <- names(typical)[sapply(typical, is.numeric)]
    if (missing(statistics)) {
        typical <- object$output[[1]][[1]]  ## first scenario, first replicate
        stat1 <- stat0
        stat2 <- c('RB','RSE','COV')
    }
    else {
        stat1 <- statistics[statistics %in% stat0]
        stat2 <- statistics[statistics %in% c('true','RB','RSE','COV','ERR')]
    }
    if (any(stat2 %in% c('true','RB','RSE','COV','ERR')) & (estname == '')) {
        stat2 <- character(0)
        warning ("cannot compute requested statistics with your data")
    }

    extractfn <- function (out, true, estimated) {
        getij <- function(df, i, j) {
            if (nrow(df) == 0)
                rep(NA, length(j))
            else
                df[i,j]
        }
        if (length(stat1)>0) {
            tmp <- lapply(out, getij, parameter, stat1)
            tmp <- do.call (rbind, tmp)
            rownames(tmp) <- 1:nrow(tmp)
        }
        else
            tmp <- matrix(nrow = length(out), ncol = 0)

        for (st in stat2) {
            if (st == 'true') {
                tmp <- cbind(tmp, rep(true, nrow(tmp)))
            }
            if (st == 'RB') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    tmp <- cbind(tmp, (est - true) / true)
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
            if (st == 'RSE') {
                est <- sapply(out, getij, parameter, estname)
                SE.est <- sapply(out, getij, parameter, SEname)
                tmp <- cbind(tmp, SE.est/est)
            }
            if (st == 'ERR') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    tmp <- cbind(tmp, abs(est-true))
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
            if (st == 'COV') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    lcl <- sapply(out, getij, parameter, 'lcl')
                    ucl <- sapply(out, getij, parameter, 'ucl')
                    tmp <- cbind(tmp, as.numeric((true>lcl) & (true<ucl)))
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
        }
        colnames(tmp) <- c(stat1, stat2)
        tmp
    }
    estimated <- sapply(object$scenarios$fitindex,
                        function(x) {
                            method <- if ('method' %in% names(object$fit.args))
                                object$fit.args$method
                            else
                                object$fit.args[[x]]$method
                            no <- object$fit & (method != 'none')
                            ifelse(length(no) == 0, TRUE, no)
                        }
                        )
    trueD <- object$scenarios[,'D']    ## vector length = number of scenarios

    if (object$outputtype == 'regionN') {
        true <- trueD   ## object$scenarios[,'D']
        true <- true * attr(object, 'regionarea')  ## for each scenario
    }
    else if (object$outputtype %in% c('predicted','derived')){
        if (parameter == 'D')
            true <- trueD
        else
            true <- object$scenarios[,parameter]
    }
    else true <- NA
    object$output <- mapply(extractfn,
                            object$output,
                            true,
                            estimated,
                            SIMPLIFY = FALSE)
    object$outputtype <- 'numeric'
    class(object) <- c('selectedstatistics', 'secrdesign', 'list')
    attr(object, 'parameter') <- parameter
    object
}

make.array <- function (object) {
    inputs <- attr(object$scenarios, 'inputs')
    if (is.null(inputs))
        stop("array output requires 'inputs' attribute")
    dims <- sapply(inputs, length)
    vdims <- dims[dims>1]
    varying <- names(inputs)[dims > 1]
    statistics <- dimnames(object$output[[1]])[[2]]
    nstat <- length(statistics)
    nrepl <- nrow(object$output[[1]])
    if (length(varying)>0) {
        outputarray <- array (dim = c(nrepl, nstat, vdims), dimnames =
                              c(list(1:nrepl), list(statistics), inputs[varying]))
    }
    else {
        outputarray <- array (dim = c(nrepl, nstat), dimnames =
                              list(NULL, statistics))
    }
    outputarray[] <- unlist(object$output)
    outputarray
}

predict.fittedmodels <- function (object, ...) {
    output <- lapply(object$output, lapply, predict, ...)
    object$output <- output
    object$outputtype <- 'predicted'
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

coef.fittedmodels <- function (object, ...) {
    output <- lapply(object$output, lapply, coef, ...)
    object$output <- output
    object$outputtype <- 'coef'
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

derived.SL <- function (object, ...) {
    if (!inherits(object,'fittedmodels'))
        stop ("require fitted secr models")
    output <- lapply(object$output, lapply, derived, ...)
    object$output <- output
    object$outputtype <- 'derived'
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

regionN.SL <- function (object, ...) {
    if (!inherits(object,'fittedmodels'))
        stop ("require fitted secr models")
    output <- lapply(object$output, lapply, region.N, ...)
    ra <- sapply(output, function(x) attr(x[[1]], 'regionarea'))
    object$output <- output
    object$outputtype <- 'regionN'
    attr(object, 'regionarea') <- ra ## one per scenario
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}


argdf <- function (args) {
    if (is.null(args))
        NULL
    else {
        tmp <- lapply(args, function(x) sapply(x, format))
        nm <- unique(unlist(lapply(tmp, names)))
        tmp0 <- matrix('', nrow = length(tmp), ncol = length(nm),
                       dimnames = list(NULL,  nm))
        for (i in 1:length(tmp)) tmp0[i,names(tmp[[i]])] <- tmp[[i]]
        as.data.frame(tmp0)
    }
}

header <- function (object) {
    values <- lapply(object$scenarios, unique)
    nvalues <- sapply(values, length)
    notID <- names(object$scenarios) != 'scenario'
    constant <- object$scenarios[1, (nvalues==1) & notID]
    rownames(constant) <- 'value'
    constant <- data.frame(t(constant))
    constant[,1] <- as.character(constant[,1])
    varying <- object$scenarios[, nvalues>1, drop = FALSE]

    ti <- unique(object$scenarios$trapsindex)
    detectors <- data.frame(trapsindex = ti, trapsname = names(object$trapset)[ti])

    pi <- unique(object$scenarios$popindex)
    popargs <- argdf(object$pop.args)
    if (!is.null(popargs))
        popargs <- data.frame(popindex = pi, popargs[pi,,drop = FALSE])

    di <- unique(object$scenarios$detindex)
    detargs <- argdf(object$det.args)
    if (!is.null(detargs))
        detargs <- data.frame(detindex = di, detargs[di,,drop = FALSE])

    fi <- unique(object$scenarios$fitindex)
    fitargs <- argdf(object$fit.args)
    if (!is.null(fitargs))
        fitargs <- data.frame(fitindex = fi, fitargs[fi,,drop = FALSE])

    list(call = object$call, starttime = object$starttime, proctime = object$proctime,
         constant = constant, varying = varying, pop.args = popargs, det.args = detargs,
         fit.args = fitargs, detectors = detectors, nrepl = object$nrepl,
         outputclass = class(object))
}

summary.secrdesign <- function (object, ...) {
    out <- list(header = header(object))
    class (out) <- c('summarysecrdesign', 'list')
    out
 }

summary.rawdata <- function (object, ...) {
    tmp <- fit.models(object)
    summary(tmp, ...)
 }

summary.estimatetables <- function (object, ...) {
    tmp <- select.stats(object)
    summary(tmp, ...)
}

summary.selectedstatistics <- function (object, dec = 5,
                                fields = c('n', 'mean', 'se'),
                                alpha = 0.05, type = c('list','dataframe','array'),
                                ...) {
    if (length(fields) == 1)
        if (tolower(fields) == 'all')
            fields <- c('n', 'mean', 'sd', 'se', 'min', 'max', 'lcl', 'ucl', 'q025',
                        'median', 'q975', 'rms')

    z      <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    alpha1 <- 0.025
    alpha3 <- 0.975
    qfields <- fields[substring(fields,1,1)=='q']
    if (length(qfields)>0) alpha1 <- as.numeric(substring(qfields[1],2,4))/1000
    if (length(qfields)>1) alpha3 <- as.numeric(substring(qfields[2],2,4))/1000
    q1tag <- paste('q', formatC(1000*alpha1, width=3, format="d", flag = "0"), sep='')
    q3tag <- paste('q', formatC(1000*alpha3, width=3, format="d", flag = "0"), sep='')

    type <- match.arg(type)

    ## allow for zero-length
    minNA <- function (x) if (sum(!is.na(x)) > 0) min(x, na.rm = T) else NA
    maxNA <- function (x) if (sum(!is.na(x)) > 0) max(x, na.rm = T) else NA

    ## data.frame output
    if (type %in% c('list','dataframe')) {
        sumx <- function (x) {
            n    <- sum(!is.na(x))
            mean <- mean(x, na.rm = T)
            sd   <- sd(x, na.rm = T)
            se   <- sd/n^0.5
            minx <- minNA(x)
            maxx <- maxNA(x)
            lcl  <- mean - z * se
            ucl  <- mean + z * se
            q    <- quantile(x, na.rm = T, prob = c(alpha1, 0.5, alpha3))
            rms  <- mean(x^2, na.rm = T)^0.5
            tmp  <- c(n, mean, sd, se, minx, maxx, lcl, ucl, q, rms)
            tmp[is.nan(tmp)] <- NA
            names(tmp) <- c('n', 'mean', 'sd', 'se', 'min', 'max', 'lcl', 'ucl', q1tag,
                            'median', q3tag, 'rms')
            tmp
        }
        valstring <- function (slist) {
            slist <- t(slist)   ## to reorder fields
            val <- as.numeric(slist)
            names(val) <- apply(do.call(expand.grid, dimnames(slist)), 1, paste,
                                collapse='.')
            val
        }
        scenariolabel <- function (scen) {
            sv <- scen[varying]
            sv <- sapply(sv, format)
            nv <- names(object$scenarios)[varying]
            label <- paste(nv, sv, sep=' = ')
            paste(label, collapse=', ')
        }
        tidy <- function(x) {
            x <- as.data.frame(x)
            x[,fields, drop=F]
        }
        varying <- apply(object$scenarios,2, function(x) length(unique(x)))>1
        varying[1] <- FALSE ## drop scenario code
        tmp <- lapply (object$output, function(y) round(t(apply(y,2,sumx)), dec))
        tmp <- lapply(tmp, tidy)

        if (type == 'list') {
            names(tmp) <- apply(object$scenario,1,scenariolabel)
        }
        else {
            new <- do.call(rbind, lapply(tmp, valstring))
            tmp <- cbind(object$scenarios[,varying], new)
        }
        out <- list(header = header(object), OUTPUT = tmp)
    }
    ## array output
    else {
        outputarray <- make.array(object)
        nd     <- length(dim(outputarray))
        n      <- apply(outputarray, 2:nd, function(y) sum(!is.na(y)))
        mean   <- apply(outputarray, 2:nd, mean, na.rm = T)
        rms    <- apply(outputarray, 2:nd, function (x) mean(x^2, na.rm = T)^0.5)
        sd     <- apply(outputarray, 2:nd, sd, na.rm = T)
        se     <- sd/n^0.5
        minx   <- apply(outputarray, 2:nd, minNA)
        maxx   <- apply(outputarray, 2:nd, maxNA)
        q      <- apply(outputarray, 2:nd, quantile, na.rm = T, prob =
                        c(alpha1, 0.5, alpha3))
        out    <- list(n = n,
                       mean = mean,
                       se = se,
                       sd = sd,
                       min = minx,
                       max = maxx,
                       lcl = mean - z * se,
                       ucl = mean + z * se,
                       rms = rms,
                       q1 = apply(q,2:nd,'[',1),
                       median = apply(q,2:nd,'[',2),
                       q3 = apply(q,2:nd,'[',3))
        names(out)[names(out)=='q1'] <- q1tag
        names(out)[names(out)=='q3'] <- q3tag
        qi <- grepl('q', fields)
        if (any(qi)) {
            fields <- fields[!qi]
            fields <- c(fields, q1tag, q3tag)
        }
        out    <- out[fields]
        statlast <- function (y) {
            y[is.nan(y)] <- NA
            y <- round(y, dec)
            aperm(y, c(2:(nd-1),1))
        }
        tmp    <- lapply(out, statlast)
        out <- list(header = header(object), OUTPUT = tmp)
    }
    class (out) <- c('summarysecrdesign', 'list')
    out
}

print.summarysecrdesign <- function (x, ...) {
    print(x$header$call)
    cat('\n')
    cat('Replicates   ', x$header$nrepl, '\n')
    cat('Started      ', x$header$starttime, '\n')
    cat('Run time     ', round(x$header$proctime/60,3), ' minutes \n')
    cat('Output class ', x$header$outputclass[1], '\n')
    cat('\n')
    print(x$header['constant'])
    cat("$varying\n")
    print(x$header[['varying']], row.names = FALSE)
    cat("\n")
    cat("$detectors\n")
    print(x$header[['detectors']], row.names = FALSE)
    cat("\n")
    if (!is.null(x$header[['pop.args']])) {
        cat("$pop.args\n")
        print (x$header[['pop.args']], row.names = FALSE)
        cat("\n")
    }
    if (!is.null(x$header[['det.args']])) {
        cat("$det.args\n")
        print (x$header[['det.args']], row.names = FALSE)
        cat("\n")
    }
    if (!is.null(x$header[['fit.args']])) {
        cat("$fit.args\n")
        print (x$header[['fit.args']], row.names = FALSE)
        cat("\n")
    }
    if (!is.null(x$OUTPUT)) {
        cat('OUTPUT\n')
        print(x$OUTPUT)
    }
}

#####################################
## plot method for secrdesign object

plot.selectedstatistics <- function (x, scenarios, statistic, type = c('hist','CI'),
                                     refline, xlab = NULL, ...) {
    plothist <- function (mat, stat, scen, ...) {
        if (is.character(stat) & !(stat %in% colnames(mat)))
            stop ("requested statistic not in data")
        else if (is.numeric(stat)) {
            if (stat > dim(mat)[2])
                stop ("statistic exceeds dimension of data")
            stat <- colnames(mat)[stat]
        }
        if (is.null(xlab)) {
            if (!is.null(param))
                xlab <- paste(stat, ' (', param, ')', sep = '')
            else
                xlab <- stat
        }
        hist (mat[,stat], main = "", xlab = xlab, ...)
        if (refline) {
            if (param %in% c('E.N','R.N')) {
                true <- x$scenarios[scen, 'D']
                true <- true * attr(x, 'regionarea')[scen]
            }
            else {
                true <- x$scenarios[scen, param]
            }
            abline(v = true, col = 'red')
        }
        mtext(side=3, line = 1, paste('Scenario', scen, collapse=' '), cex = par()$cex*1.1)
    }
    plotCI <- function (mat, scen, ...) {
        estname <- 'estimate'
        if (x$outputtype == 'coef')
            estname <- 'beta'
        if (!all(c(estname,'lcl','ucl') %in% colnames(mat)))
            stop ("requires statistics 'estimate','lcl','ucl'")
        if (is.null(xlab))
            xlab <- 'Replicate'
        if (!is.null(param))
            ylab <- paste(estname, ' (', param, ')', sep='')
        nrepl <- nrow(mat)
        plot (c(1,nrepl), range(mat[,c('lcl','ucl')], na.rm = TRUE), type = 'n',
              xlab = xlab, ylab = ylab, ...)
        segments (1:nrepl, mat[,'lcl'], 1:nrepl, mat[,'ucl'])
        if (refline) {
            if (param %in% c('E.N','R.N')) {
                true <- x$scenarios[scen, 'D']
                true <- true * attr(x, 'regionarea')[scen]
            }
            else {
                true <- x$scenarios[scen, param]
            }
            abline(h = true, col = 'red')
        }
        points (1:nrepl, mat[,estname], pch = 21, bg = 'white')
        mtext(side=3, line = 1, paste('Scenario', scen, collapse=' '), cex = par()$cex*1.1)
    }
    param <- attr(x, 'parameter')
    type <- match.arg(type)
    if (missing(scenarios))
        scenarios <- 1:length(x$output)
    if (missing(statistic))
        statistic <- 1
    if (missing(refline))
        refline <- type == 'CI'   ## default TRUE if CI, FALSE otherwise
    for (scenario in scenarios) {
        if (type == 'CI')
            plotCI(x$output[[scenario]], scen = scenario,...)
        else {
            for (s in statistic) {
                if (type == 'hist')
                    plothist(x$output[[scenario]], stat = s, scen = scenario,...)
            }
        }
    }
}
