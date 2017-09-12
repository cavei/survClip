\name{pathwaySurvivalTest}
\alias{pathwaySurvivalTest}
\title{Survival analysis on the whole pathways}
\description{This function allows to perform servival analysis on pathways. The analysis can be either topological or classical. The analysis is based on data reduction based on Principal Component Analysis and cox survival on the most influent PCs.}
\usage{
pathwaySurvivalTest(expr, survAnnot, graph,
pcsSurvCoxMethod = c("regular", "topological", "sparse"),
alwaysShrink = FALSE, maxPCs = 10, survFormula = "Surv(days, status) ~")
}
\arguments{
  \item{expr}{expression matrix.}
  \item{survAnnot}{a data frame for survival annotations specified accordin to the formuala. The data frame must contain days and status.}
  \item{graph}{a graphNEL object for a graph}
  \item{pcsSurvCoxMethod}{a method to perform pca. A method to perform "regular", "topological", "sparse" for regular pca, topological based pca and sparse pca, respectively. The latter one (sparse) is particularly suited for cliques only.}
  \item{alwaysShrink}{if TRUE, alway shrink the covariance matrix. Deafult=FALSE}
  \item{maxPCs}{maximum number of PCs to contemplate in the cox formula "Surv(days, status) ~ pc1.."}
  \item{survFormula}{the formula to use in coxoh analysis. Defaut="Surv(days, status) ~". Please note that the formula end with '~' meaning that PCs will be added. Alternatively the formula should end with '+' because the PCs with be always pasted.}
}
\details{
  With regular method, a regular pca analysis is used to compute PCs. With topological covariance matrix is penalyze according to the topology of the pathways. The max number of PCs to consider is estimated by "estim_ncp" and a maximum number of PCs can be fixed by the user.
}
\value{A survPath object.}

\author{
Paolo Martini
}
\seealso{
\code{\link{cliqueSurvivalTest}}
}
