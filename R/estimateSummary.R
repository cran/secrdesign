###############################################################################
## package 'secrdesign'
## estimateSummary.R
## 2023-02-25
## 2023-03-10 abind::adrop when no group in any scenario
###############################################################################

# Convert output (nested list of estimate tables) to array
# with dimensions (Parameter, statistic, Group, Scenario, Replicate)
estimateArray <- function (object) {
    nscen <- length(object$output)
    nrepl <- length(object$output[[1]])
    groups <- unique(object$scenarios$group)
    ngrp  <- max(1, length(groups))
    typical <- object$output[[1]][[1]]
    if (inherits(typical, "secrlist")) stop("estimateArray expects estimate table(s), not secrlist")
    if(length(dim(typical))!=2) {
        groupnames <- if (length(typical) == ngrp) groups else names(typical)
        ngrp <- length(groupnames)
        typical <- typical[[1]]
    }
    else {
        ngrp <- 1
        groupnames <- 1
        npar <- unique(unlist(lapply(object$output, sapply, nrow)))
        if (length(npar)>1) stop("number of parameters varies among scenarios or replicates")
    }
    typical <- data.frame(typical)  # convert matrix output from modelAverage
    # count actual tables per replicate of each scenario
    ntab <- sapply(object$output, 
        function(x) {
            if (length(dim(x[[1]]))==2) 1 else length(x[[1]])
        })
    if (length(table(ntab)) > 1) 
        stop("varying number of estimatetables per scenario")
    if (ntab[1] > ngrp)
        warning("multiple tables per group")
    num <- sapply(typical, is.numeric)  # copes with matrix
    if (inherits(object, 'estimatetables')) {
        values <- unlist(object$output)
        values <- suppressWarnings(as.numeric(values))  # link text -> NA
        arr <- array(values, 
            dim = c(dim(typical), ngrp, nrepl, nscen),
            dimnames = list(
                Parameter = rownames(typical),
                names(typical),
                Group = groupnames,
                Replicate = 1:nrepl,
                Scenario = 1:nscen))
        arr <- arr[,num,,, ,drop = FALSE]   # drop link field
        arr <- aperm(arr, c(1,2,3,5,4))
        arr
        
    }
    else stop("unsuitable data; expecting estimatetables from secrdesign")
}
###############################################################################

# Direct summary of output from run.scenarios
estimateSummary <- function (object, parameter = 'D', 
    statistics = c('true', 'nvalid', 'EST', 'seEST', 'RB', 
        'seRB', 'RSE', 'RMSE', 'rRMSE', 'COV'), true, validrange = c(0, Inf), 
    checkfields = c('estimate','SE.estimate'),
    format = c('data.frame', 'list'), cols = NULL) {
    arr <- estimateArray(object)
    if (missing(true)) {
        # assumes groups, if present, are used in model
        # otherwise incompatible
        true <- if (parameter %in% names(object$scenarios))
            object$scenarios[,parameter] 
        else rep(NA, nrow(object$scenarios))  # 2023-04-15 for robustness
    }
    if (length(true) != prod(dim(arr)[3:4])) stop ("incompatible 'true' vector")
    # check for out-of-range estimates
    validate <- function(x) {
        x1 <- as.numeric(x[checkfields])
        if (!all(x1 >= validrange[1] & x1 <= validrange[2])) {
            x[] <- NA
        }
        x
    }
    if (!all(checkfields %in% dimnames(arr)[[2]]))
        warning("checkfields not found in object$output; check not performed")
    else
        arr <- aperm(apply(arr, c(1,3,4,5), validate), c(2,1,3,4,5))

    # working subsets of data
    estcol <- c('estimate','SE.estimate')
    # switch for coef()
    if (all(c('beta','SE.beta') %in% dimnames(arr)[[2]])) {
        estcol <- c('beta','SE.beta') 
    }
    parm <- arr[parameter,estcol[1],,,, drop = FALSE]
    parmse <- arr[parameter,estcol,,,, drop = FALSE]
    parmcl <- arr[parameter, c('lcl','ucl'),,,, drop = FALSE]
    
    se <- function (x) sd(x, na.rm = TRUE) / sum(!is.na(x))^0.5
    nvalid <- apply(parm, 3:4, function(x) sum(!is.na(x)))
    EST <- apply(parm, 3:4, mean, na.rm = TRUE)
    seEST <- apply(parm, 3:4, se)
    
    truemat <- array(true, dim=dim(arr)[3:4], dimnames = dimnames(arr)[3:4])
    # relative bias
    rb <- function (x, true) (x-true)/true
    RB <- sweep(parm, MARGIN = 3:4, STATS = true, FUN = rb)
    mnRB <- apply(RB, 3:4, mean, na.rm = TRUE)
    seRB <- apply(RB, 3:4, se)
    
    # relative SE
    rse <- function (x) x[2]/x[1]
    RSE <- apply(parmse, 3:5, rse)
    RSE <- apply(RSE, 1:2, mean, na.rm = TRUE)
    
    # relative RMSE
    err <- function (x, true) (x-true)
    rms <- function (x) sqrt(mean(x^2, na.rm = TRUE))
    ERR <- sweep(parm, MARGIN = 3:4, STATS = true, FUN = err)
    RMSE <- apply(ERR, 3:4, rms)
    rRMSE <- RMSE / true
    
    # coverage of CI
    inrange <- function(x) (x[1,'lcl',,,drop=FALSE] < true) & (x[1,'ucl',,,drop=FALSE] > true)
    OK <- apply(parmcl, 5, inrange)
    dim(OK) <- dim(parmcl)[3:5]
    COV <- apply(OK, 1:2, mean, na.rm = TRUE)
    COV <- array(COV, dim = dim(parmcl)[3:4], dimnames = dimnames(parm)[3:4])
    
    # assemble output, discarding unwanted summaries
    out <- list(
        true = truemat, 
        nvalid = nvalid, 
        EST = EST, 
        seEST = seEST, 
        RB = mnRB, 
        seRB = seRB, 
        RSE = RSE, 
        RMSE = RMSE, 
        rRMSE = rRMSE, 
        COV = COV)[statistics]
    if (names(out)[1] == 'true') names(out)[1] <- paste0("true.", parameter)
    
    ngrp <- max(1, length(unique(object$scenarios$group)))
    if (dim(arr)[3] == 1) {
        out <- lapply(out, abind::adrop, 1)
        if (is.null(cols)) {
            cols <- 'scenario'
        }
    }
    else {
        if (is.null(cols)) {
            if ('group' %in% names(object$scenarios) && ngrp > 1)
                cols <- c('scenario', 'group')
            else
                cols <- 'scenario'
        }
    }
        
    # optionally recast output list as data.frame
    format <- match.arg(format)
    if (format == 'data.frame') {
        out  <- do.call(data.frame, lapply(out, as.numeric))
        # pre-pend columns from object$scenarios
        if (dim(arr)[3]>1 && ngrp>1) {
            # valid groups defined in scenarios df
            rows <- 1:nrow(object$scenarios)   
        }
        else {
            # groups from multi-session, etc.
            if (dim(arr)[3]>1) {
                rows <- rep(unique(object$scenarios$scenario), each = dim(arr)[3])
            }
            else  {
                rows <- 1:nrow(out)
            }
            rows <- match(rows, object$scenarios$scenario)
        }
        scen <- object$scenarios[rows,cols, drop = FALSE]
        if (!is.null(cols)) out <- cbind(scen, out)
    }
    
    out
}
