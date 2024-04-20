# test dynamic traps option
# 2022-10-19

library(secrdesign)

sigma <- 25
box <- data.frame(x = c(0,0,1,1,0), y = c(0,1,1,0,0))
box[] <- box * 10 * sigma

scen <- make.scenarios(D = 5, sigma = sigma, lambda0 = 0.2, detectfn = 14)
traplist <- list(trapfn = make.lacework)
trapargs <- list(list(region = box, spacing = c(100, 20), origin = NULL))
masklist <- list(make.mask(box, buffer = 5*sigma, nx = 32))

sims <- run.scenarios(nrepl = 10, scenarios = scen, trapset = traplist, 
    maskset = masklist, trap.args = trapargs, ncores = 1, seed = 1237)

test_that("Random-origin lacework counts as expected", {
    expect_equal(summary(sims)$OUTPUT[[1]]$mean, 
        c(31.60000, 18.70000, 15.40000, 1.47398, 0.23584, 21.03075), tolerance = 1e-4)
})
