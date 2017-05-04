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

test_survivalClip <- function(){
  set.seed(1234)
  test <- survClip(expr, survAnnot, graph)
}

