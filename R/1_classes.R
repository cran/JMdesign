setClass("powerLongSurv",
    representation(
        title = "character",
        subtitle = "character",
        t= "vector",
        p= "vector",
        N="integer",
        nevents="integer",
        censr="numeric",
        tmedian="numeric",
        meantf="numeric",
        SigmaTheta="matrix",
        ordtraj="integer",
        BSigma="matrix",
        beta = "numeric",
        alpha="numeric",
        power="numeric"
    )
)
