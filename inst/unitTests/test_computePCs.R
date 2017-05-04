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
isto <- ann[,"Isto"]
grade <- ann[,"Grado"]

survAnnot <- data.frame(days=time, status=events, grade=grade,
                        isto=isto, row.names=samples)

k <- pathways("hsapiens", "kegg")
k <- convertIdentifiers(k, "entrez")
graph <- pathwayGraph(k[["Pathways in cancer"]])

genes <- intersect(graph::nodes(graph), row.names(exp))
graph <- graph::subGraph(genes, graph)
expr <- exp[genes, , drop=FALSE]

test_computePCsTopo <- function(){
  cliques <- clipper:::extractCliquesFromDag(graph, root=NULL)
  test <- computePCs(t(expr), npc=1, robust=FALSE, shrink=FALSE, cliques=cliques)
  checkEqualsNumeric(dim(test), c(73,1))
}

test_computePCsTopo2Pc <- function(){
  cliques <- clipper:::extractCliquesFromDag(graph, root=NULL)
  test <- computePCs(t(expr), npc=2, robust=FALSE, shrink=FALSE, cliques=cliques)
  checkEqualsNumeric(dim(test), c(73,2))
}

test_computePCs <- function(){
  test <- computePCs(t(expr), npc=1, robust=FALSE, shrink=FALSE, cliques=NULL)
  checkEqualsNumeric(dim(test), c(73,1))
}

test_computePCsRobust <- function(){
  test <- computePCs(t(expr), npc=1, robust=TRUE, shrink=FALSE, cliques=NULL)
  checkEqualsNumeric(dim(test), c(73,1))
}
