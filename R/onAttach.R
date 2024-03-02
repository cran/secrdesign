###############################################################################
## package 'secr'
## onAttach.R
## last changed 2024-03-02
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('secrdesign'), ' ', .local$packageType)
    packageStartupMessage( "This is secrdesign ", version, 
                           ". For overview type ?secrdesign" )
    options(matprod = "internal")  # in case OpenBLAS is a problem
}
