#' Class \code{"powerLongSurv"}
#' 
#' Class of objects like the output of function \code{"powerLongSurv()"}.
#'
#' @docType class
#'
#' @section Objects from the Class:
#'
#'  Objects can be created by calls of the form \code{new("powerLongSurv", ...)}.
#'
#' @slot title Object of class \code{"character"}
#' @slot subtitle Object of class \code{"character"}
#' @slot t Object of class \code{"vector"}
#' @slot p Object of class \code{"vector"}
#' @slot N Object of class \code{"integer"}
#' @slot nevents Object of class \code{"integer"}
#' @slot censr Object of class \code{"numeric"}
#' @slot tmedian Object of class \code{"numeric"}
#' @slot meantf Object of class \code{"numeric"}
#' @slot SigmaTheta Object of class \code{"matrix"}
#' @slot ordtraj Object of class \code{"integer"}
#' @slot BSigma Object of class \code{"matrix"}
#' @slot beta Object of class \code{"numeric"}
#' @slot alpha Object of class \code{"numeric"}
#' @slot power Object of class \code{"numeric"}
#' 
#' @method show powerLongSurv
#'
#' @author Emil A. Cornea, Liddy M. Chen, Bahjat F. Qaqish, Haitao Chu, and 
#'   Joseph G. Ibrahim
#'
#' @seealso \code{\link{powerLongSurv}}, \code{\link{show-methods}}
#'
#' @examples 
#' showClass("powerLongSurv")
#'
#' @import methods
#' @keywords classes
#' @export
setClass("powerLongSurv",
    representation(
        title = "character",
        subtitle = "character",
        t = "vector",
        p = "vector",
        N = "integer",
        nevents = "integer",
        censr = "numeric",
        tmedian = "numeric",
        meantf = "numeric",
        SigmaTheta = "matrix",
        ordtraj = "integer",
        BSigma = "matrix",
        beta = "numeric",
        alpha = "numeric",
        power = "numeric"),
    prototype(
      title = "No Calculation Performed",
      subtitle = "No Calcultion Performed",
      t = numeric(1),
      p = numeric(1),
      N = 0L,
      nevents = 0L,
      censr = 0.0,
      tmedian = 0.0,
      meantf = 0.0,
      SigmaTheta = matrix(0,1,1),
      ordtraj = 0L,
      BSigma = matrix(0,1,1),
      beta = 0.0,
      alpha = 0.0,
      power = 0.0)
)
