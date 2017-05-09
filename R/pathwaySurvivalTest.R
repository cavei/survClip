pathwaySurvivalTest <- function(expr, survAnnot, graph, nperm=100, perc=0.8, formula="Surv(days, status)~pc", root=NULL, alwaysShrink){
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
  shrink <- length(samples) < maxcliques || alwaysShrink

  days   <- survAnnot$days
  events <- survAnnot$status

  gtpvalue           <- globaltest::p.value(globaltest::gt(Surv(days, events==1), alternative=t(expr), permutations=nperm))
  pcspvalue          <- survCoxOnAllPCs(genes, expr, perc, survAnnot, formula, pc2class=TRUE)
  pcspvalueCov       <- survCoxOnAllPCs(genes, expr, perc, survAnnot, formula, pc2class=TRUE, shrink=FALSE, cliques=cliques)
  pcspvalueCovAlways <- survCoxOnAllPCs(genes, expr, perc, survAnnot, formula, pc2class=TRUE, shrink=TRUE, cliques=cliques)

  return(list(gtPvalue=gtpvalue, pcsPvalue=pcspvalue, pcsPvalueCov=pcspvalueCov, pcsPvalueCovAlways=pcspvalueCovAlways))
}
