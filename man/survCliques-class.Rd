\name{survCliques-class}
\docType{class}
\alias{survCliques-class}
\alias{survCliques}

\title{Class \code{"survCliques"}}
\description{Class for survival clique analysis}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("survCliques", ...)}.
alphas, zlist, cliques, cliquesLoadings.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}}
    \item{\code{alphas}:}{Object of class \code{"numeric"}. It contains
      pvalues of the cliques.}
    \item{\code{zlist}:}{Object of class \code{"list"}. For each
      cliques, the list of pvalues for all the covariates.}
    \item{\code{cliques}:}{Object of class \code{"list"}. The list of
      the genes in the cliques.}
    \item{\code{coxObjs}:}{Object of class \code{"list"}. For each
      cliques, the object used for Coxph analysis.}
    \item{\code{cliquesLoadings}:}{Object of class \code{"list"}. For each
      cliques, the loadings as calculated by the PCA.}
    \item{\code{cliquesExpr}:}{Object of class \code{"list"}. For each
      cliques, the cliques expression.}
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
