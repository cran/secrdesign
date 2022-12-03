## 2022-10-21 start

library(secrdesign)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################
trm <- make.grid(8, 8, detector = "multi")
trp <- make.grid(8, 8, detector = "proximity")
msk <- make.mask(trm, buffer = 100, type = 'trapbuffer')

nrm.m <- Enrm(D = 5, trm, msk, list(lambda0 = 0.2, sigma = 20), 5)
nrm.p <- Enrm(D = 5, trp, msk, list(lambda0 = 0.2, sigma = 20), 5)

# nrm.m
# En       Er       Em 
# 20.42885 32.59759 28.62447 
# nrm.p
# En       Er       Em 
# 20.42885 56.14689 49.64859 

Qpm.p <- Qpm(D = 5, trp, msk, list(lambda0 = 0.2, sigma = 20), 5)

en2.m <- En2(D = 5, trm, msk, list(lambda0 = 0.2, sigma = 20), 5)
en2.p <- En2(D = 5, trp, msk, list(lambda0 = 0.2, sigma = 20), 5)

test_that("Enrm correct for 'multi'", {
    expect_equal(unname(nrm.m), c(20.42885, 32.59759, 28.62447), tolerance = 1e-4)
})

test_that("Enrm correct for 'proximity'", {
    expect_equal(unname(nrm.p), c(20.42885, 56.14689, 49.64859), tolerance = 1e-4)
})

test_that("Qpm correct for 'proximity'", {
    expect_equal(unname(Qpm.p), c(0.3817752, 0.2816757), tolerance = 1e-4)
})

test_that("En2 correct for 'multi'", {
    expect_equal(unname(en2.m), c(20.42885, 13.80840), tolerance = 1e-4)
})

test_that("En2 correct for 'proximity'", {
    expect_equal(unname(en2.p), c(20.42885, 15.07251), tolerance = 1e-4)
})



###############################################################################
