library(graphite)
library(survival)
data(exp)
data(ann)

set.seed(1234)

samples <- intersect(colnames(exp),row.names(ann))

exp <- exp[,samples, drop=FALSE]
ann <- ann[samples,, drop=FALSE]

time <- computeDays(ann[,1:2])
events <- as.numeric(ann[,3])

survAnnot <- data.frame(days=time, status=events, row.names=samples)

k <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(k[["Pathways in cancer"]], "entrez")
graph <- pathwayGraph(p)

genes <- intersect(graph::nodes(graph), row.names(exp))
graph <- graph::subGraph(genes, graph)
expr <- exp[genes, , drop=FALSE]

test_computePCs_topological <- function(){
  cliques <- clipper:::extractCliquesFromDag(graph, root=NULL)
  test <- computePCs(t(expr), shrink=FALSE, method="topological", cliques=cliques, maxPCs=10)
  checkEqualsNumeric(dim(test$x), c(73,10)) ##
}

test_computePCs_regular <- function(){
  test <- computePCs(t(expr), shrink=FALSE, method="regular", maxPCs=10)
  checkEqualsNumeric(dim(test$x), c(73,10)) ##
}

test_computePCs_sparse <- function(){
  test <- computePCs(t(expr), shrink=FALSE, method="sparse", maxPCs=10)
  checkEqualsNumeric(dim(test$x), c(73,10)) ##
}
