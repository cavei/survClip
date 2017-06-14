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

test_cliqueSurvivalTest <- function(){
	set.seed(1234)
	test <- cliqueSurvivalTest(exp, graph, survAnnot)
	checkTrue(test@alphas[10] <= 0.05)
}
