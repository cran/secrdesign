##############################################################################
## package 'secrdesign'
## transformOutput.R
## 2023-04-15
##############################################################################

transformOutput <- function (object, extractfn, outputtype = "predicted", ...) {
    object$output <- lapply(object$output, function(x) lapply(x, extractfn, ...))
    outputtype(object) <- outputtype
    object$extractfn <- extractfn
    object
}
