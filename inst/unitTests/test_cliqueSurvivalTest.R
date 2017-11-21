library(graphite)
library(survival)
data(exp)
data(survAnnot)
data(graph)

row.names(exp) <- paste0("ENTREZID:", row.names(exp))

test_cliqueSurvivalTest <- function(){
	set.seed(1234)
	test <- cliqueSurvivalTest(exp, survAnnot, graph)
	checkTrue(test@alphas[10] <= 0.05)
}
