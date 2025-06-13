library(secrdesign)

# simulate some capture data
set.seed(124)
grid <- make.grid(8, 8, spacing = 2, detector = 'count')
pop <- sim.popn(D = 1000, core = grid, buffer = 4)
ch <- simOU.capthist(grid, pop, detectpar = list(tau = 50, sigma = 1, epsilon = 0.25), 
                     noccasions = 100, savepath = TRUE, by = 10)

test_that("Counts correct for 'simOU.capthist'", {
    expect_equal(nrow(ch), 23, tolerance = 1e-4)
    expect_equal(sum(ch), 120, tolerance = 1e-4)
})

test_that("RPSV correct for 'simOU.capthist'", {
    expect_equal(RPSV(ch, CC = TRUE), 0.2973264, tolerance = 1e-6)
})


