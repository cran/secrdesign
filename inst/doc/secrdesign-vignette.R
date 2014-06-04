### R code from vignette source 'secrdesign-vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: secrdesign-vignette.Rnw:42-46
###################################################
options(continue=" ")
options(width=70)
opar <- par()
on.exit (par(opar))


###################################################
### code chunk number 2: secrdesign-vignette.Rnw:135-140 (eval = FALSE)
###################################################
## library(secrdesign)
## scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
## traps1 <- make.grid()
## sims1 <- run.scenarios(nrepl = 50, trapset = traps1, scenarios =
##      scen1, seed = 345, fit = TRUE)


###################################################
### code chunk number 3: secrdesign-vignette.Rnw:142-144
###################################################
library(secrdesign)
load('sims1.RData')


###################################################
### code chunk number 4: secrdesign-vignette.Rnw:151-152
###################################################
summary(sims1)$OUTPUT


###################################################
### code chunk number 5: secrdesign-vignette.Rnw:181-185
###################################################
library(secrdesign)
mydetectors <- list(grid6x6 = make.grid(6,6),
                    grid8x9 = make.grid(8,9),
                    grid12x12 = make.grid(12,12))


###################################################
### code chunk number 6: secrdesign-vignette.Rnw:213-215
###################################################
make.scenarios (trapsindex = 1:3, noccasions = 4, D = 5, g0 = 0.2,
                sigma = c(20,30))


###################################################
### code chunk number 7: secrdesign-vignette.Rnw:231-233
###################################################
make.scenarios (trapsindex = 1:3, noccasions = c(8,4,2), D = 5, g0 = 0.2,
                sigma = c(20,30), crosstraps = FALSE)


###################################################
### code chunk number 8: secrdesign-vignette.Rnw:399-401 (eval = FALSE)
###################################################
## run.scenarios (nrepl = 1000, scenarios = scen1, trapset = traps1,
##     extractfn = closedN, estimator = c("null", "h2", "jackknife"))


###################################################
### code chunk number 9: secrdesign-vignette.Rnw:583-586
###################################################
stats1 <- select.stats(sims1, parameter = "D", statistics = c("estimate",
   "lcl", "ucl", "RB", "RSE", "COV"))
lapply(stats1$output, head, 4)


###################################################
### code chunk number 10: secrdesign-vignette.Rnw:706-707
###################################################
summary(stats1)


###################################################
### code chunk number 11: secrdesign-vignette.Rnw:722-725 (eval = FALSE)
###################################################
## par(mfrow = c(2,2))
## plot(stats1, type = "hist", statistic = "estimate")
## plot(stats1, type = "CI")


###################################################
### code chunk number 12: secrdesign-vignette.Rnw:783-786
###################################################
sims2 <- run.scenarios(nrepl = 50, trapset = traps1, scenarios = scen1,
    fit = TRUE, fit.args = list(method = "none"))
summary(sims2)


###################################################
### code chunk number 13: secrdesign-vignette.Rnw:842-851 (eval = FALSE)
###################################################
## scen3 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
## traps3 <- make.grid()
## raw3 <- run.scenarios(nrepl = 20, trapset = traps3, scenarios =
##     scen3, fit = FALSE, extractfn = identity)
## summary(raw3)
## ## fit and summarise models
## sims3 <- fit.models(raw3, fit.args = list(list(model = g0~1),
##     list(model = g0~T)), fit = TRUE, ncores = 4)
## summary(sims3)


