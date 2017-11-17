library(graphite)
library(survival)
data(exp)

row.names(exp) <- paste0("ENTREZID:", row.names(exp))
k <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(k[["Pathways in cancer"]], "entrez")
graph <- pathwayGraph(p)

genes <- intersect(graph::nodes(graph), row.names(exp))
graph <- graph::subGraph(genes, graph)
expr <- exp[genes, , drop=FALSE]

test_computePCs_topological <- function(){
  cliques <- clipper:::extractCliquesFromDag(graph, root=NULL)
  test <- computePCs(t(expr), shrink=FALSE, method="topological", cliques=cliques, maxPCs=5)
  checkEqualsNumeric(dim(test$x), c(73,5)) ##
}

test_computePCs_regular <- function(){
  test <- computePCs(t(expr), shrink=FALSE, method="regular", maxPCs=5)
  checkEqualsNumeric(dim(test$x), c(73,5)) ##
}

test_computePCs_sparse <- function(){
  test <- computePCs(t(expr), shrink=FALSE, method="sparse", maxPCs=5)
  checkEqualsNumeric(dim(test$x), c(73,5)) ##
}
