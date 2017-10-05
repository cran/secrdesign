##############################################################################
## package 'secrdesign'
## getdetectpar.R
## 2017-07-21 moved from Lambda.R
##############################################################################

getdetectpar <- function (D, C, sigma = NULL, k = 0.5, ...) {
    args <- list(...)
    args$D <- D
    if (is.null(sigma)) {
        if (!is.null(args$detectfn))
            if (secr:::valid.detectfn(args$detectfn, valid = 14:18) != 14)
                stop ("must specify sigma for detectfn not 'HHN'")
        sigma <- 100 * k / sqrt(D)
    }
    if (missing(C)) stop("C (total number of detections)  must be provided")
    objective <- function (lambda0) {
        if (is.null(args$detectpar))
            args$detectpar <- list(sigma=sigma, lambda0 = lambda0)
        else {
            args$detectpar$sigma <- sigma
            args$detectpar$lambda0 <- lambda0
        }
        C - sum(do.call(Enrm, args)) 
    }
    L0 <- uniroot(objective, interval=c(1e-4, 1e4))
    if (L0$estim.prec > 0.001) warning ("error in lambda0 exceeds 0.001")
    out <- list(lambda0 = L0$root, sigma = sigma)
    if (!is.null(args$detectfn)) {  ## otherwise default to HHN
        if (args$detectfn == "HHR") out$z <- args$detectpar$z
        if (args$detectfn == "HAN") out$w <- args$detectpar$w
    }
    out
}

