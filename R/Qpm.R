##############################################################################
## package 'secrdesign'
## Qpm.R
## 2022-10-23, 2022-11-30
##############################################################################

Qpm <- function (D, traps, mask, detectpar, noccasions, detectfn = 
        c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX')) {
    if (min(dist(traps)) == 0) 
        warning("not all detector locations unique")
    if (!(detector(traps)[1] %in% c('proximity','count')))
        stop ("Qpm is for detector types proximity and count only")
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:18))
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    detectpars <- unlist(detectpar[secr:::parnames(detectfn)])
    D <- rep(D, length.out = nrow(mask)) * attr(mask, 'area')  # per cell
    temp <- Qpmcpp (
        as.double(unlist(detectpars)), 
        as.double(D),
        as.matrix(edist(traps, mask)), 
        as.integer(detectfn),
        as.integer(noccasions)
    )
    
    c(Qp = temp$Qp, Qpm = temp$Qpm)
    
}
