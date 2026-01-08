##############################################################################
## package 'secrdesign'
## En2.R
## 2022-10-23, 2022-11-19, 2022-11-30
##############################################################################

En2 <- function (D, traps, mask, detectpar, noccasions, detectfn = 
        c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX')) {
    if (min(dist(traps)) == 0) 
        warning("not all detector locations unique")
    if (!(detector(traps)[1] %in% c('multi','proximity','count')))
        stop ("En2 is for detector types multi, proximity and count only")
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::secr_valid.detectfn(detectfn, valid = c(0,1,2,14:19))
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    detectpars <- unlist(detectpar[secr:::secr_parnames(detectfn)])
    dettype <- secr:::secr_detectorcode(traps, noccasions = noccasions)[1]
    D <- rep(D, length.out = nrow(mask)) * attr(mask, 'area')  # per cell
    temp <- En2cpp (
        as.integer(dettype),
        as.double(unlist(detectpars)), 
        as.double(D),
        as.matrix(edist(traps, mask)), 
        as.integer(detectfn),
        as.integer(noccasions)
    )
    
    c(En = temp$En, En2 = temp$En2)
    
}
