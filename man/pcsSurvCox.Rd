\name{pcsSurvCox}
\alias{pcsSurvCox}

\title{Cox Analysis on Principal Components}
\description{Survival Analysis using Cox Model on Principal Components. For internal use only.}
\usage{
pcsSurvCox(genes, expr, annotations,
method = c("regular", "topological", "sparse"),
shrink = FALSE, cliques = NULL, maxPCs = 10,
survFormula = "Surv(days, status) ~")
}
\arguments{
  \item{genes}{vector of genes}
  \item{expr}{whole expression matrix}
  \item{annotations}{survival annotations}
  \item{method}{PCA methods}
  \item{shrink}{shrink covariance matrix}
  \item{cliques}{when method=topological use this cliques to define topology in IPF}
  \item{maxPCs}{maximum number of PC to consider}
  \item{survFormula}{survival forumula to use}
}
\details{
for internal use only
}
\value{cox model object}
\author{Paolo Martini}
