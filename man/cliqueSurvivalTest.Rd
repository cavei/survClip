\name{cliqueSurvivalTest}
\alias{cliqueSurvivalTest}
\title{Perform survival test on all the cliques of the graph.}
\description{This function performs survival test on given pathway using
  a matrix with survival annotation.}
\usage{cliqueSurvivalTest(expr, survAnnot, graph, pcsSurvCoxMethod = c("regular", "sparse"),
alwaysShrink = FALSE, maxPCs = 10, survFormula = "Surv(days, status) ~")}
\arguments{
  \item{expr}{expression matrix}
  \item{graph}{a graphNEL object}
  \item{survAnnot}{a data frame for survival annotations specified according to the survFormuala. The data frame must contain days and status.}
  \item{pcsSurvCoxMethod}{a method to perform PCA. A method to perform
    "regular", "sparse" for regular PCA and sparse PCA, respectively. }
  \item{alwaysShrink}{if TRUE, always shrink the covariance matrix. Deafult=FALSE}
  \item{maxPCs}{maximum number of PCs used in the Cox formula "Surv(days, status) ~ PC1.."}
  \item{survFormula}{the formula to use in coxph analysis. Defaut="Surv(days, status) ~". Please note that the formula end with '~' meaning that PCs will be added.}
}
\details{Survival test is made according to survFormula. With 'regular' method, a regular PCA analysis is used to compute
  PCs. With 'sparse' method, a penalized regression is used for the
  estimation of PCs (as implemented in elasticnet).}
\value{A survCliques object.}
\author{Paolo Martini}
\seealso{
  \code{\link{pathwaySurvivalTest}}
  \code{\link{getTopLoadGenes}}
}

\examples{
if (require(graphite)) {
  data(exp)
  data(survAnnot)
  data(graph)
  row.names(exp) <- paste0("ENTREZID:", row.names(exp))
  genes <- intersect(graph::nodes(graph), row.names(exp))
  graph <- graph::subGraph(genes, graph)
  expr <- exp[genes, , drop=FALSE]
  cliqueSurvivalTest(expr, survAnnot, graph, maxPCs=2)
}
}