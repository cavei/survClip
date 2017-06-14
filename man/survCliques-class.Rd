\name{survCliques-class}
\docType{class}
\alias{survCliques-class}
\alias{survCliques}

\title{Class \code{"survCliques"}}
\description{Descriptor for survival clique analysis}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("survCliques", ...)}.
alphas=alphas, zlist=zlist, cliques=cliques
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}}
    \item{\code{alphas}:}{Object of class \code{"numeric"} pvalues }
    \item{\code{zlist}:}{Object of class \code{"list"} zlist }
    \item{\code{cliques}:}{Object of class \code{"list"} cliques }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{Ulist}"}, by class "list", distance 2.
Class \code{"\linkS4class{Uvector}"}, by class "list", distance 3.
}
\section{Methods}{
No methods defined with class "survCliques" in the signature.
}
\author{
Paolo Martini
}
\examples{
showClass("survCliques")
}
\keyword{classes}
