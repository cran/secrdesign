## ----setup, eval = TRUE, echo=FALSE, message=FALSE---------------------------------
options(width = 85)
makepng <- FALSE

## ----Ex1out, eval = TRUE-----------------------------------------------------------
summary(sims1)$OUTPUT

## ----layouts, eval = TRUE----------------------------------------------------------
library(secrdesign)
mydetectors <- list(grid6x6 = make.grid(6,6),
                    grid8x9 = make.grid(8,9),
                    grid12x12 = make.grid(12,12))

## ----make.s, eval = FALSE----------------------------------------------------------
#  make.scenarios (trapsindex = 1, noccasions = 3, nrepeats = 1, D, g0,
#      sigma, lambda0, detectfn = 0, recapfactor = 1, popindex = 1,
#      detindex = 1, fitindex = 1, groups, crosstraps = TRUE)

## ----scen1, eval = TRUE------------------------------------------------------------
make.scenarios (trapsindex = 1:3, noccasions = 4, D = 5, g0 = 0.2, sigma = c(20,30))

## ----scen2, eval = TRUE------------------------------------------------------------
make.scenarios (trapsindex = 1:3, noccasions = c(8,4,2), D = 5, g0 = 0.2,
                sigma = c(20,30), crosstraps = FALSE)

## ----run.sc, eval = FALSE----------------------------------------------------------
#  run.scenarios (nrepl, scenarios, trapset, maskset, xsigma = 4,
#      nx = 32, pop.args, det.args, fit = FALSE, fit.args, chatnsim = 0,
#      extractfn = NULL, multisession = FALSE, ncores = 1, byscenario = TRUE,
#      seed = 123)

## ----owncod, eval = FALSE----------------------------------------------------------
#  sum1 <- function(out) {
#      require(abind)
#      ## collapse replicates to an array, omitting non-numeric column
#      out <- do.call(abind, c(out, along = 3))[,-1,,drop = FALSE]
#      ## convert array from character to numeric
#      mode(out) <- "numeric"
#      ## take the average over replicates (meaningless for some fields)
#      apply(out, 1:2, mean, na.rm = TRUE)
#      }
#  lapply(closedNsim$output, sum1)

## ----selsta, eval=FALSE------------------------------------------------------------
#  select.stats(object, parameter = "D", statistics)

## ----fndpar, eval=FALSE------------------------------------------------------------
#  find.param(object)

## ----seldem, eval = TRUE-----------------------------------------------------------
stats1 <- select.stats(sims1, parameter = "D", statistics = c("estimate",
   "lcl", "ucl", "RB", "RSE", "COV"))
lapply(stats1$output, head, 4)

## ----valid, eval = FALSE-----------------------------------------------------------
#  x <- validate (x, test, validrange = c(0, Inf), targets = test)

## ----summry, eval = FALSE----------------------------------------------------------
#  summary (object, dec = 5, fields = c("n", "mean", "se"), alpha = 0.05,
#      type = c("list", "dataframe", "array"), ...)

## ----stafld,echo=FALSE, comment=''-------------------------------------------------
bullet <-  '*' ## rawToChar(as.raw(149))
tmp <- matrix(bullet, nr=8, ncol=12)
dimnames(tmp) <- list(Statistics = c("estimate","SE.estimate","lcl","ucl","RB","RSE","ERR","COV"),
Fields = c("n","mean","se","sd","min","max","lcl","ucl","rms","median","q025","q975"))
tmp[1:6,9] <- ""
tmp[3:5,7:8] <- ""
tmp[8,7:12] <- ""
print(tmp, quote = FALSE)

## ----Ex1sum, eval = TRUE-----------------------------------------------------------
summary(stats1, c('n', 'mean', 'se', 'median'))

## ----plot, eval=FALSE--------------------------------------------------------------
#  par(mfrow = c(2,2))
#  plot(stats1, type = "hist", statistic = "estimate")
#  plot(stats1, type = "CI")

## ----plotpng, echo=FALSE, eval=makepng---------------------------------------------
#  png(file='d:/density secr 3.1/secrdesign/vignettes/secrdesign-fig3.png',
#      width = 850, height = 800)
#  par(mfrow = c(2,2), cex=1.2)
#  plot(stats1, type = "hist", statistic = "estimate")
#  plot(stats1, type = "CI")
#  dev.off()

## ----Ex2Sum, eval = TRUE-----------------------------------------------------------
summary(sims2)

## ----extrfn------------------------------------------------------------------------
defaultextractfn <- function(x) {
    if (inherits(x, 'try-error')) {
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
    }
    else if (inherits(x, 'capthist')) {
        ## summarised raw data
        counts <- function(CH) {
            ## for single-session CH
            if (nrow(CH)==0) { ## 2015-01-24
                if (sighting(traps(CH))) 
                    c(n = 0, ndet = 0, nmov = 0, dpa = 0,
                      unmarked=0, nonID = 0, nzero = 0) 
                else
                    c(n=0, ndet=0, nmov=0, dpa = NA)
            }
            else {
                n <- nrow(CH)
                ndet <- sum(abs(CH)>0)
                r <- ndet - n
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
                else
                    c(n=n, ndet=ndet, nmov=nmoves, dpa = dpa, rse = sqrt(1/n + 1/r))
            }
        }
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
    else if (inherits(x,'secr') & (!is.null(x$fit))) {
        ## fitted model:
        ## default predictions of 'real' parameters
        out <- predict(x)
        if (is.data.frame(out))
            out
        else {
            warning ("summarising only first session, group or mixture class")
            out[[1]]
        }
    }
    else
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
}

## ----Ex4Plt, eval = FALSE, fig.width = 8, fig.height = 8---------------------------
#  par(mfrow = c(4,3), cex=0.8, mgp=c(2.2,0.6,0))
#  plot(sims4, statistic = "n", breaks = seq(0,80,5))      ## animals
#  plot(sims4, statistic = "nmov", breaks = seq(0,140,5))  ## movements

## ----Ex4png, eval=makepng, echo=FALSE----------------------------------------------
#  png(file='d:/density secr 3.1/secrdesign/vignettes/secrdesign-fig4.png',
#      width = 850, height = 800)
#  par(mfrow = c(4,3), cex=1)
#  plot(sims4, statistic = "n", breaks = seq(0,80,5))  ## number of animals
#  plot(sims4, statistic = "nmov", breaks = seq(0,140,5))
#  dev.off()

## ----Ex5Sum------------------------------------------------------------------------
## select statistics and throw out any replicates with SE > 100, if any
stats5 <- select.stats(sims5)
stats5 <- validate(stats5, "SE.estimate", c(0,100), "all")
sum5 <- summary(stats5, fields = c("n","mean","se","lcl","ucl", "median"))

## ----Ex5Plt------------------------------------------------------------------------
## plot
plot(c(0.5,6.5), c(-0.2,0.4), type = "n", xlab = "Scenario", ylab = "RB(D-hat)")
for (i in 1:12) {
    xv <- if (i<=6) i else (i-6)+0.05
    segments (xv, sum5$OUTPUT[[i]]["RB","lcl"], xv, sum5$OUTPUT[[i]]["RB","ucl"])
    ptcol <- if (i<=6) "white" else "black"
    points(xv, sum5$OUTPUT[[i]]["RB","mean"], pch = 21, bg = ptcol)
}
abline(h = 0, col="red")
text(c(1.5,3.5,5.5), rep(0.38,3), paste("recapfactor", c(0.5,1,2), sep = " = "))

## ----Ex5Out, eval = TRUE-----------------------------------------------------------
## look at extended output
sum5

## ----Ex6Sum------------------------------------------------------------------------
summary(sims6)

## ----Ex6a--------------------------------------------------------------------------
sims6a <- run.scenarios (1, scen6, traps(possumCH), possummask,
    pop.args = poplist, det.args = list(savepopn = TRUE),
    extractfn = identity)

## ----Ex6apl, eval = TRUE, fig.width=8.5, fig.height = 4----------------------------
## sims6a$output is now a list (one component per scenario) of lists
## (one component per replicate) of simulated capthist objects, each
## with its 'popn' object embedded as an attribute

pop1 <- attr(sims6a$output[[1]][[1]], "popn")
pop2 <- attr(sims6a$output[[2]][[1]], "popn")
par(mfrow = c(1,2), mar=c(1,1,1,6))
plot(possummask, covariate = "D1", dots = FALSE, breaks = 0:6)
plot(traps(possumCH), detpar = list(col = 'green', pch = 15), add = TRUE)
plot(pop1, frame = FALSE, add = TRUE, col = "blue", pch = 16, cex = 0.6)
plot(possummask, covariate = 'D2', dots = FALSE, breaks = 0:6)
plot(traps(possumCH), detpar = list(col = 'green', pch = 15), add = TRUE)
plot(pop2, frame = FALSE, add = TRUE, col = "blue", pch = 16, cex = 0.6)

## ----Ex6aIn, eval = FALSE----------------------------------------------------------
#  ## click on map to display height; Esc to exit
#  spotHeight(possummask, prefix = "D2")

## ----Ex7, warning = FALSE, message = FALSE-----------------------------------------
library(secrlinear)

## ----Ex7Sum------------------------------------------------------------------------
summary(sims7)

