###############################################################################
## package 'secrdesign'
## saturation.R
## 2017-11-02
###############################################################################

saturation <- function (traps, mask, detectpar, detectfn = 
                        c('HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HN', 'HR', 'EX'), 
                        D, plt = FALSE, add = FALSE, ...) {
    if (!(detector(traps)[1] %in% c('multi','proximity', 'capped', 'count')))
        stop ("only for 'multi','proximity', 'capped' or 'count' detectors")
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::secr_valid.detectfn(detectfn, valid = c(0,1,2,14:19))
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 1 to 15, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    
    cellarea <- attr(mask, 'area') 
    dk <- edist(traps, mask)                                 ## K x M
    if (detectfn == 14)
        h <- exp(-dk^2/2/detectpar$sigma^2)                  ## K x M
    else if (detectfn == 15)
        h <- 1 - exp(-(dk/detectpar$sigma)^-detectpar$z)
    else if (detectfn == 16)
        h <- exp(-dk^2/2/detectpar$sigma^2)  
    else if (detectfn == 17)
        h = exp(- (dk-detectpar$w)^2 / 2 / detectpar$sigma^2)
    else if (detectfn == 18)
        h = pgamma(dk, shape = detectpar$z, scale = detectpar$sigma/detectpar$z, 
                   lower.tail = FALSE, log.p = FALSE)
    h <- detectpar$lambda0 * h
    if (detector(traps)[1] == "multi") {
        Hi <- apply(h, 2, sum)                                ## M hazard of detn | x 
        hmult <- (1 - exp(-Hi)) / Hi
        pkx <- sweep(h, MARGIN = 2, STATS = hmult, FUN = "*") ## K x M Pr caught @ k | x
        hkx <- -log(1-pkx)
    }
    else {
        hkx <- h
    }
    if (length(D) > 1)
        hkx <- sweep(hkx, MARGIN = 2, STATS = D, FUN = "*") ## K x M 
    else
        hkx <- hkx * D
    Hk <- apply(hkx, 1, sum) * cellarea       ## K  
    p <- 1-exp(-Hk)  
    out <- list(bydetector = p, mean = mean(p)) # lambda0bias = mean(p)/mean(Hk)-1
    if (plt) {
        covariates(traps) <- data.frame(saturation = 1-exp(-Hk))
        plot(as.mask(traps), covariate = 'saturation', dots = FALSE, add = add, ...)
        plot(traps, add = TRUE)
        invisible(out)
    }
    else {
        out
    }
}

