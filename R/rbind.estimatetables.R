#############################################################################
## package 'secrdesign'
## rbind.estimatetables.R
## 2023-01-12 work in progress
## 2023-08-30 c.estimatetables revised
#############################################################################

rbind.estimatetables <- function (..., deparse.level = 1) {
    # combine 2 or more objects output from run.scenarios
    
    allargs <- list(...)
    inputnames <- as.character(match.call(expand.dots=FALSE)$...) 
    
    scenarios <- lapply(allargs, '[[', 'scenarios')
    headers <- lapply(allargs, header)
    names1 <- names(scenarios[[1]])
    estnames <- lapply(allargs, function(x)  names(x$output[[1]][[1]]))
    
    nruns <- length(scenarios)
    if (nruns>1)
        for (i in 2:nruns) {
            if (any (names(scenarios[[i]] != names1))) stop ("differing columns in scenarios") 
            if (any (scenarios[[i]] != scenarios[[1]])) stop ("differing scenarios") 
            if (any (headers[[i]]$fit.args != headers[[1]]$fit.args)) stop ("differing fit.args")
            if (any (estnames[[i]] != estnames[[1]])) stop ("differing estimatetable names")
            # also check trapset, nx, det.args?
        }
    
    output <- lapply(allargs, '[[', 'output') 
    temp <- allargs[[1]]
    nscenarios <- nrow(temp$scenarios)
    for (i in seq_len(nscenarios)) {
        temp$output[[i]] <- do.call(c, lapply(output, '[[', i))
    }
    
    temp$proctime <- sum(sapply(allargs, '[[', 'proctime'))
    temp$nrepl <- sum(sapply(allargs, '[[', 'nrepl'))
    temp
}
############################################################################################
rbind.selectedstatistics <- function (..., deparse.level = 1) {
    # combine 2 or more objects output from run.scenarios
    allargs <- list(...)
    # inputnames <- as.character(match.call(expand.dots=FALSE)$...) 
    
    scenarios <- lapply(allargs, '[[', 'scenarios')
    headers <- lapply(allargs, header)
    names1 <- names(scenarios[[1]])
    estnames <- lapply(allargs, function(x)  names(x$output[[1]][[1]]))
    
    nruns <- length(scenarios)
    if (nruns>1)
        for (i in 2:nruns) {
            if (any (names(scenarios[[i]] != names1))) stop ("differing columns in scenarios") 
            if (any (scenarios[[i]] != scenarios[[1]])) stop ("differing scenarios") 
            if (any (headers[[i]]$fit.args != headers[[1]]$fit.args)) stop ("differing fit.args")
            if (any (estnames[[i]] != estnames[[1]])) stop ("differing estimatetable names")
            # also check trapset, nx, det.args?
        }
    
    output <- lapply(allargs, '[[', 'output') 
    temp <- allargs[[1]]
    nscenarios <- nrow(temp$scenarios)
    for (i in seq_len(nscenarios)) {
        temp$output[[i]] <- do.call(rbind, c(lapply(output, '[[', i)))
    }
    
    temp$proctime <- sum(sapply(allargs, '[[', 'proctime'))
    temp$nrepl <- sum(sapply(allargs, '[[', 'nrepl'))
    temp
}
############################################################################################

c.estimatetables <- function (...) {
    # combine 2 or more objects output from run.scenarios
    
    allargs <- list(...)

    headers   <- lapply(allargs, header)
    scenarios <- lapply(allargs, '[[', 'scenarios')
    output    <- lapply(allargs, '[[', 'output') 
    estnames  <- lapply(allargs, function(x)  names(x$output[[1]][[1]]))

    # check scenario compatibility
    inputscen <- sapply(scenarios, nrow) 
    names1 <- names(scenarios[[1]])
    nruns <- length(scenarios)
    if (nruns>1) {
        for (i in 2:nruns) {
            if (any (names(scenarios[[i]] != names1))) stop ("differing columns in scenarios") 
        }
    }
    
    # construct output
    temp <- allargs[[1]]   # inherit header from first object
    temp$scenarios <- do.call(rbind, scenarios)
    temp$output    <- do.call(c, output)
    temp$proctime  <- sum(sapply(allargs, '[[', 'proctime'))
    temp$nrepl     <- sapply(allargs, '[[', 'nrepl') # separate, not sum 2023-08-30
    temp
}
###############################################################################

c.selectedstatistics <- function (...) {
    
    # combine 2 or more objects output from run.scenarios
    
    allargs <- list(...)
    
    headers   <- lapply(allargs, header)
    scenarios <- lapply(allargs, '[[', 'scenarios')
    output    <- lapply(allargs, '[[', 'output') 
    estnames  <- lapply(allargs, function(x)  names(x$output[[1]][[1]]))
    
    # check scenario compatibility
    inputscen <- sapply(scenarios, nrow) 
    names1 <- names(scenarios[[1]])
    nruns <- length(scenarios)
    if (nruns>1) {
        for (i in 2:nruns) {
            if (any (names(scenarios[[i]] != names1))) stop ("differing columns in scenarios") 
        }
    }
    
    # construct output
    temp <- allargs[[1]]
    temp$scenarios <- do.call(rbind, scenarios)
    temp$output    <- do.call(c, output)
    temp$proctime  <- sum(sapply(allargs, '[[', 'proctime'))
    temp$nrepl     <- sapply(allargs, '[[', 'nrepl')
    temp
}
###############################################################################

# library(secrdesign)
# setwd('d:/density communication/ipsecr/single catch simulations')
# sims1 <- readRDS(paste0('originalsims', 'ML', 1, '.RDS' ))
# sims2 <- readRDS(paste0('originalsims', 'ML', 2, '.RDS' ))
# source('d:/density secr 4.5/secrdesign/R/rbind.estimatetables.R')
# sims <- rbind(sims1, sims1)
# summary(sims)
# 
# sims <- c(sims1, sims2)
# summary(sims)
