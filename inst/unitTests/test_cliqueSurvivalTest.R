library(graphite)
library(survival)
data(exp)
data(survAnnot)

k <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(k[["Pathways in cancer"]], "entrez")
graph <- pathwayGraph(p)

test_cliqueSurvivalTest <- function(){
	set.seed(1234)
	test <- cliqueSurvivalTest(exp, graph, survAnnot)
	checkTrue(test@alphas[10] <= 0.05)
}
