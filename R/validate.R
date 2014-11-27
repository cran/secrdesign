## 2014-04-17, 27, 2014-10-30
validate <- function (x, test, validrange = c(0, Inf), targets = test) {
    onescenario <- function (scen, minx, maxx) {
        criterion <- (scen[,test] < minx) | (scen[,test] > maxx)
        criterion[is.na(criterion)] <- TRUE     ## 2014-10-30
        scen[criterion,targets] <- NA
        scen[is.na(scen[,test]),targets] <- NA
        scen
    }
    if (!inherits(x, 'selectedstatistics'))
        stop ("requires 'selectedstatistics'")
    if (tolower(targets) == "all")
        targets <- colnames(x$output[[1]])
    if (!all(c(test, targets) %in% colnames(x$output[[1]])))
        stop ("test or targets not found in input")
    nscen <- length(x$output)
    if (!is.matrix(validrange))
        validrange <- matrix (validrange, byrow = TRUE, nrow = nscen, ncol = 2)
    if (nrow(validrange) != nscen)
        stop ("invalid `validrange'")
    minx <- validrange[,1]
    maxx <- validrange[,2]
    x$output <- mapply(onescenario, x$output, minx, maxx, SIMPLIFY = FALSE)
    x
}


