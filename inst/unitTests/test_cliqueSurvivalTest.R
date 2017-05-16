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

test_cliqueSurvivalTest <- function(){
	set.seed(1234)
	test <- cliqueSurvivalTest(exp, graph, survAnnot)
	checkTrue(test$alpha[10] <= 0.05)
}






