## ----setup, eval = TRUE, echo=FALSE, message=FALSE---------------------------------
library(secrdesign)
options(width = 85)

## ----costarg, echo = FALSE---------------------------------------------------------
str(costing)

## ----pervis1-----------------------------------------------------------------------
covariates(tr) <- data.frame(costpervisit = 5 + edist(tr, tr[1,])/20)

## ----arglist-----------------------------------------------------------------------
str(optimalSpacing)

## ----osout-------------------------------------------------------------------------
out$rotRSE[-1]

## ----osout2------------------------------------------------------------------------
out2$simRSE[-1]  # do not show each simulation

## ----plottrials, warning = FALSE, fig.width = 8, fig.height = 4.5------------------
plottrials <- function (trials, title = '', hline = c(0.1, 0.2, 0.5, 1), 
                        xmax, type = 'r', xl, label = '') {
    plot(0,0,type = 'n',xlim = c(0,xmax), ylim = c(0.05,2), log = 'y', 
         yaxs = 'i', xaxs = 'i',
         ylab = expression(paste('RSE(', hat(italic(D)), ')')),
         xlab = '', axes = FALSE)
    axis(1)
    axis(2, at = c(0.05,0.1,0.2,0.5,1,2), las = 1,
         labels = c('0.05','0.1','0.2','0.5','1','2'))
    box()
    lines(1:xmax, 1/sqrt(1:xmax))
    abline(h = hline, lty = 2, xpd = FALSE)
    tmp <- as.data.frame(do.call(rbind,trials))
    mtext(side = 1, line = 2.4, xl, cex = 0.9, xpd = TRUE)
    xval <- if (type == 'r') tmp$r else pmin(tmp$n,tmp$r)
    points(xval, tmp$CV, pch=21, bg = "red", cex=1.4, xpd=TRUE)
    usr <- par()$usr 
    text ((usr[2]-usr[1])*-0.2, (usr[4]-usr[3])*1.9, label, cex = 1.3, xpd = TRUE) 
    invisible()
}

par(mfrow = c(1,2), mar = c(4,5,3,1), mgp = c(2.4,0.8,0))
plottrials(trials2, '', xmax = 270, xl = expression(italic(r)), type = 'r', label = 'a.')
plottrials(trials2, '', xmax = 120, xl = expression(paste("min(",italic(n), "," , 
    italic(r),")")), type = 'nr', label = 'b.')

