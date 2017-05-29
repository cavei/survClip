pathwaySurvivalTest <- function(expr, survAnnot, graph, nperm=100, perc=0.8, formula="Surv(days, status)~pc",
                                root=NULL, alwaysShrink=FALSE, maxPCs=10){
  genes <- nodes(graph)
  genes <- intersect(genes, rownames(expr))
  if (length(genes) <= 3){
    return(NULL)
    warning("Too few genes.")
  }

  samples <- intersect(colnames(expr),row.names(survAnnot))
  if (length(sample) == 0){
    return(NULL)
    warning("No sample intersection.")
  }

  survAnnot <- survAnnot[samples,]
  expr <- expr[genes,samples, drop=FALSE]

  graph <- graph::subGraph(genes, graph)
  expr <- expr[genes,, drop=FALSE]

  cliques <- clipper:::extractCliquesFromDag(graph, root=root)

  maxcliques <- max(sapply(cliques, length))
  shrink <- length(samples) < maxcliques | alwaysShrink

  days   <- survAnnot$days
  events <- survAnnot$status

  set.seed(1234)
  gtpvalue           <- globaltest::p.value(globaltest::gt(Surv(days, events==1), alternative=t(expr), permutations=nperm))
  set.seed(1234)
  pcspvalue          <- survCoxOnAllPCs(genes, expr, perc, survAnnot, pc2class=TRUE, robust=FALSE, maxPCs=maxPCs)
  set.seed(1234)
  pcspvalueCov       <- survCoxOnAllPCs(genes, expr, perc, survAnnot, pc2class=TRUE, robust=FALSE, shrink=FALSE, cliques=cliques, maxPCs=maxPCs)
  set.seed(1234)
  pcspvalueCovAlways <- survCoxOnAllPCs(genes, expr, perc, survAnnot, pc2class=TRUE, robust=FALSE, shrink=TRUE,  cliques=cliques, maxPCs=maxPCs)
  set.seed(1234)
  pcspvalueCovAlwaysNT <- survCoxOnAllPCs(genes, expr, perc, survAnnot, pc2class=TRUE, robust=FALSE, shrink=FALSE, cliques=cliques, maxPCs=maxPCs, useTopology = FALSE)
  set.seed(1234)
  pcspvalueCovAlwaysNT.shrink <- survCoxOnAllPCs(genes, expr, perc, survAnnot, pc2class=TRUE, robust=FALSE, shrink=TRUE,  cliques=cliques, maxPCs=maxPCs, useTopology = FALSE)
  
  return(list(gtPvalue=gtpvalue, pcsPvalue=pcspvalue, pcsPvalueCov=pcspvalueCov, pcsPvalueCovAlways=pcspvalueCovAlways, noTopo=pcspvalueCovAlwaysNT, noTopoShrinked=pcspvalueCovAlwaysNT.shrink))
}
