\name{survivalcox}
\alias{survivalcox}
\title{cox model analysis}
\description{Cox analysis. For internal use only}
\usage{
survivalcox(coxObj, formula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coxObj}{data.frame: patients x covariates.}
  \item{formula}{formula to use}
}
\details{
for internal use only.
}
\value{
A list with
\item{pvalue}{pvalue of the model}
\item{zlist}{pvalues of single covariates}
\item{coxObj}{the original coxObj passed to the function}
}
\author{Paolo Martini}
