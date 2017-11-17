\name{cliqueSurvivalTest}
\alias{cliqueSurvivalTest}
\title{Perform survival test on all the cliques of the graph.}
\description{Given a pathway and a matrix with survival annotation compute survival test.}
\usage{cliqueSurvivalTest(expr, survAnnot, graph, pcsSurvCoxMethod = c("regular", "sparse"),
alwaysShrink = FALSE, maxPCs = 10, survFormula = "Surv(days, status) ~")}
\arguments{
  \item{expr}{expression matrix}
  \item{graph}{a graphNEL object}
  \item{survAnnot}{a data frame for survival annotations specified accordin to the formuala. The data frame must contain days and status}
  \item{pcsSurvCoxMethod}{a method to perform PCA. A method to perform "regular", "sparse" for regular PCA and sparse PCA, respectively.}
  \item{alwaysShrink}{if TRUE, always shrink the covariance matrix. Deafult=FALSE}
  \item{maxPCs}{maximum number of PCs used in the cox formula "Surv(days, status) ~ PC1.."}
  \item{survFormula}{the formula to use in coxph analysis. Defaut="Surv(days, status) ~". Please note that the formula end with '~' meaning that PCs will be added. If you want to include further covariates please note that the formula should end with '+' because the PCs will be always pasted at the end.}
}
\details{Survival test is made according to formula.}
\value{A survCliques object.}
\author{Paolo Martini}
\seealso{
\code{\link{pathwaySurvivalTest}}
}

\examples{
if (require(graphite)) {
  data(exp)
  data(survAnnot)
  row.names(exp) <- paste0("ENTREZID:", row.names(exp))
  k <- pathways("hsapiens", "kegg")
  p <- convertIdentifiers(k[["Pathways in cancer"]], "entrez")
  graph <- pathwayGraph(p)
  genes <- intersect(graph::nodes(graph), row.names(exp))
  graph <- graph::subGraph(genes, graph)
  expr <- exp[genes, , drop=FALSE]
  cliqueSurvivalTest(expr, survAnnot, graph, maxPCs=2)
}
}