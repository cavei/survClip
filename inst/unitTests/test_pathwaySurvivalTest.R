library(graphite)
library(survival)
data(exp)
data(survAnnot)
data(graph)

row.names(exp) <- paste0("ENTREZID:", row.names(exp))
genes <- intersect(graph::nodes(graph), row.names(exp))
graph <- graph::subGraph(genes, graph)
expr <- exp[genes, , drop=FALSE]

test_pathwaySurvivalTest <- function(){
  set.seed(1234)
  test <- pathwaySurvivalTest(expr, survAnnot, graph, maxPCs=2)
  checkTrue(length(test@zlist) == 2)
}

