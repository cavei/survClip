\name{pathwaySurvivalTest}
\alias{pathwaySurvivalTest}
\title{Survival analysis on the whole pathway}
\description{This function performs survival analysis on pathways. The analysis can be either topological or classical. The analysis is based on data reduction based on Principal Component Analysis and a Cox proportional hazard model on the most influent PCs.}
\usage{
pathwaySurvivalTest(expr, survAnnot, graph,
pcsSurvCoxMethod = c("regular", "topological", "sparse"),
alwaysShrink = FALSE, maxPCs = 10, survFormula = "Surv(days, status) ~")
}
\arguments{
  \item{expr}{expression matrix.}
  \item{survAnnot}{a data frame for survival annotations specified according to the survFormula. The data frame must contain days and status.}
  \item{graph}{a graphNEL object for a graph}
  \item{pcsSurvCoxMethod}{a method to perform PCA. Can be "regular", "topological", "sparse" for regular PCA, topological based PCA and sparse PCA, respectively. The latter one (sparse) is particularly suited for cliques only.}
  \item{alwaysShrink}{if TRUE, always shrink the covariance matrix. Deafult=FALSE}
  \item{maxPCs}{maximum number of PCs used in the cox formula "Surv(days, status) ~ PC1.."}
  \item{survFormula}{the formula used in Coxph analysis. Defaut="Surv(days, status) ~". Please note that the formula end with '~' meaning that PCs will be added.}
}
\details{
  With 'regular' method, a regular PCA analysis is used to compute
  PCs. With 'topological' method, the covariance matrix is estimated
  using the topology of the pathway with IPS algorithm. With 'sparse' method, a penalized regression is used for the
  estimation of PCs (as implemented in elasticnet). 
  The max number of PCs used by the model is estimated by "estim_ncp" in FactoMineR. A
  maximum number of PCs can be fixed by the user. The minimum ot the two is chosen.
}
\value{A survPath object.}

\author{
Paolo Martini
}
\seealso{
\code{\link{cliqueSurvivalTest}}
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
  pathwaySurvivalTest(expr, survAnnot, graph, maxPCs=2)
}
}
