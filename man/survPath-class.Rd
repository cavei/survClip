\name{survPath-class}
\docType{class}
\alias{survPath-class}
\alias{survPath}

\title{Class \code{"survPath"}}
\description{Descriptor for survival pathway analysis.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("survPath", ...)}.
pvalue, zlist, coxObj, loadings and the method used.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} }
    \item{\code{pvalue}:}{Object of class \code{"numeric"} that contains the pvalue of the whole models generated for the analysis }
    \item{\code{zlist}:}{Object of class \code{"numeric"} list of all singular pvalues for all the covariates}
    \item{\code{coxObj}:}{Object of class \code{"data.frame"} the object use for coxph analysis}
    \item{\code{loadings}:}{Object of class \code{"matrix"} the loadings as calculated by the PCA}
    \item{\code{method}:}{Object of class \code{"character"} type of method used. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{Ulist}"}, by class "list", distance 2.
Class \code{"\linkS4class{Uvector}"}, by class "list", distance 3.
}
\section{Methods}{
No methods defined with class "survPath" in the signature.
}
\author{
Paolo Martini
}
\examples{
showClass("survPath")
}
\keyword{classes}
