\name{powerLongSurv-class}
\Rdversion{1.1}
\docType{class}
\alias{powerLongSurv-class}
\alias{JMdesign-powerLongSurv-class}

\title{Class \code{"powerLongSurv"}}
\description{
Class of objects like the output of function \code{"powerLongSurv()"}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("powerLongSurv", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{title}:}{Object of class \code{"character"}}
    \item{\code{subtitle}:}{Object of class \code{"character"}}
    \item{\code{t}:}{Object of class \code{"vector"}}
    \item{\code{p}:}{Object of class \code{"vector"}}
    \item{\code{N}:}{Object of class \code{"integer"} }
    \item{\code{nevents}:}{Object of class \code{"integer"}}
    \item{\code{censr}:}{Object of class \code{"numeric"}}
    \item{\code{tmedian}:}{Object of class \code{"numeric"}}
    \item{\code{meantf}:}{Object of class \code{"numeric"}}
    \item{\code{SigmaTheta}:}{Object of class \code{"matrix"}}
    \item{\code{ordtraj}:}{Object of class \code{"integer"}}
    \item{\code{BSigma}:}{Object of class \code{"matrix"}}
    \item{\code{beta}:}{Object of class \code{"numeric"}}
    \item{\code{alpha}:}{Object of class \code{"numeric"}}
    \item{\code{power}:}{Object of class \code{"numeric"}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "powerLongSurv")}}
	 }
}

\author{
Emil A. Cornea, Liddy M. Chen, Bahjat F. Qaqish, Haitao Chu, and Joseph G. Ibrahim

Maintainer: Shannon T. Holloway <sthollow@ncsu.edu>
}

\seealso{
    \code{\link{powerLongSurv}}, \code{\link{show-methods}}
}
\examples{
showClass("powerLongSurv")
}
\keyword{classes}
