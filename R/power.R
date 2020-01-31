power <- function (D2D1, RSE, adjustRSE, testtype, alpha) {
    if (adjustRSE) {
        sdiff <- (log(1 + RSE^2) + log(1 + RSE^2 / D2D1))^0.5
        effect <- log(D2D1) - log(sqrt(1 + RSE^2 / D2D1)) + log(sqrt(1 + RSE^2))
    }
    else {
        sdiff <- (2 * log(1 + RSE^2))^0.5
        effect <- log(D2D1)
    }
    effect <- effect / sdiff
    if (testtype == "two.sided") {
        z.alpha <- qnorm(1-alpha/2)
        pnorm(effect - z.alpha, 0, 1, lower.tail = TRUE) +
            pnorm(-effect - z.alpha, 0, 1, lower.tail = TRUE)
    }
    else if (testtype == "decrease") {
        z.alpha <- qnorm(1-alpha)
        pnorm(-effect - z.alpha, 0, 1, lower.tail = TRUE)
    }
    else {
        z.alpha <- qnorm(1-alpha)
        pnorm(effect - z.alpha, 0, 1, lower.tail = TRUE)
    }
}

powerCI <- function (D2D1, RSE, adjustRSE, alpha) {
    ## effect D2D1 is ratio of final and initial density estimates
    ## RSE is relative standard error of initial density estimate
    ##
    ## find SD and mean of effect size on log scale
    if (adjustRSE) {
        sdiff <- (log(1 + RSE^2) + log(1 + RSE^2 / D2D1))^0.5
        effect <- log(D2D1) - log(sqrt(1 + RSE^2 / D2D1)) + log(sqrt(1 + RSE^2))
    }
    else {
        sdiff <- (2 * log(1 + RSE^2))^0.5
        effect <- log(D2D1)
    }
    ## return back-transformed Wald interval for effect on log scale
    z.alpha <- qnorm(c(alpha/2, 1-alpha/2))
    exp(effect + sdiff * z.alpha)
}

preplus <- function(x) paste0(symnum(x, c(-Inf, 0, Inf), c("", "+")), x)

plotpower <- function (RSE = seq(0.05,0.25,0.05), effectrange = c(-0.99,1.5), adjustRSE = FALSE, 
                       testtype = "two.sided", alpha = 0.05, targetpower = 0.80, 
                       effectincr = 0.02, plt = TRUE, col = topo.colors(8), lty = 1, 
                       xlab = 'Population change %', ylab = 'Power %',
                       add = FALSE, shading = 'lightgreen', ...) {
    if (!add) {
        xlim <- effectrange*100
        if (xlim[1] < -98) xlim[1] <- -100
        plot(0,0,type='n', xlim = xlim, ylim=c(0,100), yaxs='i', xaxs = 'i', axes = FALSE,
             xlab = "", ylab = "")
        at <- seq(-100, xlim[2], 50)
        lab <- preplus(at)
        axis(1, at = at, labels = lab, las = 1)
        axis(2, las = 1)
        mtext(side=1, line=2.6, xlab, xpd = TRUE)
        mtext(side=2, line=3, ylab, xpd = TRUE)
        abline(h = 100*targetpower, lty = 2, xpd = FALSE)
        box()
    }
    xval <- seq(effectrange[1], effectrange[2], effectincr)
    ## get critical values
    zero <- which.min(abs(xval))
    nRSE <- length(RSE)
    col <- rep(col,nRSE)[1:nRSE]
    lty <- rep(lty,nRSE)[1:nRSE]
    lower <- upper <- numeric(nRSE)
    dpower <- function (x, rse, target = targetpower) {
        power(x+1, rse, adjustRSE, testtype, alpha) - target
    }
    pow <- matrix(nrow = length(xval), ncol = nRSE, dimnames=list(xval, RSE))
    for (i in 1:nRSE) {
        pow[,i] <- power(xval+1, RSE[i], adjustRSE, testtype, alpha)
        powlo <- pow[1,i]
        if (!is.na(powlo) && (powlo >= targetpower)) {
            lower[i] <- uniroot(dpower, interval = xval[c(1,zero)], rse = RSE[i])$root
            lower100 <- lower[i]*100
            if (!is.na(shading)) {
                polygon (c(-100,-100,lower100, lower100), c(0,100,100, 0), col = shading)
                text (lower100, 105, as.character(round(lower100, 1)), cex=0.8, xpd = TRUE)
            }
        }
        else lower[i] <- NA
        powhi <- pow[length(xval),i]
        if (!is.na(powhi) && (powhi >= targetpower)) {
            upper[i] <- uniroot(dpower, interval = xval[c(zero, length(xval))], rse = RSE[i])$root
            upper100 <- upper[i]*100
            if (!is.na(shading)) {
                polygon (c(upper100, upper100, effectrange[2]*100, effectrange[2]*100), c(0,100,100, 0), col = shading)
                text (upper100, 105, as.character(round(upper100, 1)), cex=0.8, xpd = TRUE)
            }
        }
        else upper[i] <- NA
        
        powerpct <- 100*power(xval+1, RSE = RSE[i], adjustRSE, testtype, alpha)
        lines (xval*100, powerpct, col = col[i], lty = lty[i], ...)
    }
    
    list(RSE = RSE, effectrange = effectrange, testtype = testtype, adjustRSE = adjustRSE,
         alpha = alpha, targetpower = targetpower, lower=lower, upper = upper, power = pow)
}

plotpowerCI <- function (RSE = seq(0.05,0.25,0.05), effectrange = c(-0.99,1.5), 
                         estimatedrange = effectrange, adjustRSE = FALSE, 
                       alpha = 0.05, effectincr = 0.02, col = topo.colors(8), plt = TRUE, 
                       xlab = 'Population change %', ylab = 'Estimated population change %',
                       add = FALSE, ...) {
    
    if (!add) {
        plot(0,0, xlim = effectrange+1, ylim = estimatedrange+1, xaxs='i', yaxs='i', type='n', axes = FALSE,
             xlab = '', ylab = '')
        labx <- preplus(seq(-100, effectrange[2]*100,50))
        laby <- preplus(seq(-100, estimatedrange[2]*100,50))
        axis (1, at = seq(0,effectrange[2]+1, 0.5), labels = labx, las = 1)
        axis (2, at = seq(0,estimatedrange[2]+1, 0.5), labels = laby, las = 1)
        mtext (side=1, line=2.6, xlab)
        mtext (side=2, line=3.5, ylab)
        abline(v=1, lty=2, xpd = FALSE)
        abline(h=1, lty=2, xpd = FALSE)
        # abline(0,1, lty=2, col='blue')
        box()
    }
    xval <- seq(effectrange[1], effectrange[2], effectincr) +1
    nRSE <- length(RSE)
    ci <- array(dim=c(length(xval), 2, nRSE), dimnames = list(xval, c('lower','upper'),RSE))
    for (i in 1:nRSE) {
        ci[,,i] <- t(sapply(xval, powerCI, RSE = RSE[i], adjustRSE = adjustRSE, alpha = alpha))
        lines(xval, ci[,1,i], col = col[i], xpd = FALSE, ...)
        lines(xval, ci[,2,i], col = col[i], xpd = FALSE, ...)
    }
    
    list(RSE = RSE, effectrange = effectrange, adjustRSE = adjustRSE,
         alpha = alpha, limits = ci)
}


#  
# plotpower(effectrange=c(-0.6,1.0), adjust=T, shading=NA, lwd = 2)
# text(c(24,42,57,69,79), c(82, 71, 62, 53, 44)+2, paste0(seq(0.05,0.25,0.05)*100,'%'),
#      cex = 0.8, adj=0, col='black')
# 
# par(pty='s')
# junk <- plotpowerCI(effectrange=c(-1,1.5), adjust=T, lwd = 1.5)
# text(c(1.68,1.81,1.94,2.07,2.18), c(1.48, 1.42, 1.34, 1.25, 1.14)+0.15, paste0(seq(0.05,0.25,0.05)*100,'%'),
#      cex = 0.8, adj=0, col='black')
# 
# plotpowerCI(RSE=0.2, col='red', effectrange=c(-1,1.5), adjust=T, lwd = 1.5, add=T)


