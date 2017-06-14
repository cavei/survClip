\name{cliqueSurvivalTest}
\alias{cliqueSurvivalTest}
\title{Perform survival test on all the cliques of the graph.}
\description{Given a pathway and a matrix with annotation pData compute survival test.}
\usage{cliqueSurvivalTest(expr, graph, survAnnot, pcsSurvCoxMethod = c("regular", "sparse"), alwaysShrink = FALSE, maxPCs = 10, survFormula = "Surv(days, status) ~")}
\arguments{
  \item{expr}{expression matrix}
  \item{graph}{a graphNEL object}
  \item{survAnnot}{a data frame for survival annotations specified accordin to the formuala. The data frame must contain days and status}
  \item{pcsSurvCoxMethod}{a method to perform pca. A method to perform "regular", "sparse" for regular pca and sparse pca, respectively.}
  \item{alwaysShrink}{if TRUE, alway shrink the covariance matrix. Deafult=FALSE}
  \item{maxPCs}{maximum number of PCs to contemplate in the cox formula "Surv(days, status) ~ pc1.."}
  \item{survFormula}{the formula to use in coxoh analysis. Defaut="Surv(days, status) ~". Please note that the formula end with '~' meaning that PCs will be added. Alternatively the formula should end with '+' because the PCs with be always pasted.}
}
\details{Survival test is made according to formula.}
\value{A survCliques object.}
\author{Paolo Martini}
\seealso{
\code{\link{pathwaySurvivalTest}}
}
