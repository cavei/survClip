\name{computePCs}
\alias{computePCs}
\alias{topoCompPCs}
\alias{sparseCompPCs}
\alias{compPCs}
\title{For internal usage only. Compute principal components according to the chosen method.}
\description{For internal usage. Compute the PCs according to 3 methods: regular, ropological, sparse.}
\usage{computePCs(exp, shrink, method = c("regular", "topological", "sparse"),
cliques = NULL, maxPCs)}
\arguments{
  \item{exp}{exp}
  \item{shrink}{shrink boolean}
  \item{method}{either regular, topological, sparse}
  \item{cliques}{if NULL, no topology.}
  \item{maxPCs}{mac pcs returned.}
}
\details{For internal usage.}
\value{
x=scores, sdev=sd
\item{x}{pcs}
\item{sdev}{standar deviation for pcs.}
}
\author{Paolo Martini}