## ----plottest, results = "hold", fig.width=6, fig.height=6---------------
plotc <- function (v = 'n') {
    dtype <- c('multi','proximity','count')
    detector <- rep(dtype,2)[out$trapsindex]
    out2 <- matrix(nrow=3, ncol = 2, dimnames = list(dtype, c('mean','sd')))
    for (d in 1:3) {
        OK <- (detector == dtype[d])
        RB <- (out[OK,v] - out[OK, paste0('E',v)]) / out[OK, paste0('E', v)]
        # use log scales to spread values
        plot(out[OK, paste0('E',v)], out[OK,v], log='xy',
             xlab = 'Expected', ylab = 'simulated')
        abline(0,1) # y = x line
        # return mean and sd of estimated relative bias 
        out2[d,] <- round(c(mean = mean(RB), sd = sd(RB)),5)
        mtext (side=3, line=0.2, paste(v, " ", dtype[d]), cex = 0.9)
    }
    out2
}
par(mfrow=c(3,3), mgp=c(2.3,0.6,0), mar=c(4,4,2,1), pty='s')
cat("Relative discrepancy between expected and simulated counts\n")
cat("Number of individuals\n")
plotc('n')
cat("Number of recaptures\n")
plotc('r')
cat("Number of movements\n")
plotc('m')

