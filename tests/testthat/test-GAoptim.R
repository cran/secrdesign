## 2022-10-05 start

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

###############################################################################
test_that("random number generator stable", {
    set.seed(1235)
    expect_equal(rnorm(1), -0.6979879, tolerance = 1e-6)
})
###############################################################################

library(secrdesign)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################

## GAoptim()

msk <- make.mask(type = 'rectangular', spacing = 10, nx = 30, ny = 20, buffer = 0)
alltrpsm <- make.grid(nx = 29, ny = 19, origin = c(10,10), spacing = 10, detector = 'multi')
alltrpsp <- make.grid(nx = 29, ny = 19, origin = c(10,10), spacing = 10, detector = 'proximity')

# 10 generations for demonstration, use more in practice
opt4 <- GAoptim(msk, alltrpsm, ntraps = 20, 
    detectpar = list(lambda0 = 0.5, sigma = 20), 
    detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, 
    verbose = 0, criterion = 4, seed = 123)

opt5p <- GAoptim(msk, alltrpsp, ntraps = 20, 
    detectpar = list(lambda0 = 0.5, sigma = 20), 
    detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, 
    verbose = 0, criterion = 5, seed = 123)

test_that("Genetic algorithm on track, criterion 4, multi", {
  expect_equal(opt4$des$bestobj, -46.01154, tolerance = 1e-4)
})

test_that("Genetic algorithm on track, criterion 5, proximity", {
    expect_equal(opt5p$des$bestobj, -23.56278, tolerance = 1e-4)
})

# This test fails on Mac M1

# opt5 <- GAoptim(msk, alltrpsm, ntraps = 20,
#     detectpar = list(lambda0 = 0.5, sigma = 20),
#     detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10,
#     verbose = 0, criterion = 5, seed = 123)
# 
# test_that("Genetic algorithm on track, criterion 5, multi", {
#     expect_equal(opt5$des$bestobj, -19.01397, tolerance = 1e-4)
# })

###############################################################################
