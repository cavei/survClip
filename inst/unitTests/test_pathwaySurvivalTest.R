library(graphite)
library(survival)
data(exp)
data(survAnnot)

k <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(k[["Pathways in cancer"]], "entrez")
graph <- pathwayGraph(p)

genes <- intersect(graph::nodes(graph), row.names(exp))
graph <- graph::subGraph(genes, graph)
expr <- exp[genes, , drop=FALSE]

test_pathwaySurvivalTest <- function(){
  set.seed(1234)
  test <- pathwaySurvivalTest(expr, survAnnot, graph, maxPCs=2)
  checkTrue(length(test@pvalues) == 5)
}

