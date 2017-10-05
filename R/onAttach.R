###############################################################################
## package 'secr'
## onAttach.R
## last changed 2017-10-05
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('secrdesign'), .local$packageType)
    packageStartupMessage( "This is secrdesign ", version, 
                           ". For overview type ?secrdesign" )
}
