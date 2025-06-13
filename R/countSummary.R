countSummary <- function (object, verbose = FALSE, dec = 3) {
    output <- object$output
    if (is.null(attr(output[[1]][[1]], 'counts'))) 
        stop ("output appears not to have 'counts' attribute")
    onescenario <- function(output) t(sapply(output, attr, 'counts'))
    counts <- lapply(output, onescenario)
    if (verbose) {
        desc <- function(x) round(c(n = sum(!is.na(x)),
                              mean = mean(x, na.rm = TRUE),
                              se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                              sd = sd(x, na.rm = TRUE),
                              min = min(x, na.rm = TRUE),
                              max = max(x, na.rm = TRUE)), dec)
        lapply(counts, apply, 2, desc)
    }
    else {
        t(sapply(counts, apply, 2, mean, na.rm = TRUE))
    }
}

