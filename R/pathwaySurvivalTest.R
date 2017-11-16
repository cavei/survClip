pathwaySurvivalTest <- function(expr, survAnnot, graph, pcsSurvCoxMethod=c("regular", "topological", "sparse"),
                                alwaysShrink=FALSE, maxPCs=10, survFormula = "Surv(days, status) ~"){
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

  cliques <- clipper:::extractCliquesFromDag(graph)

  maxcliques <- max(sapply(cliques, length))
  shrink <- length(samples) < maxcliques | alwaysShrink

  days   <- survAnnot$days
  events <- survAnnot$status
  
  method = pcsSurvCoxMethod[1]
  
  res <- pcsSurvCox(genes, expr, survAnnot, method=method, shrink=alwaysShrink, cliques=cliques, maxPCs=maxPCs, survFormula = survFormula)
  new("survPath",
      pvalue = res$pvalue, zlist = res$zlist, coxObj = res$coxObj, loadings = res$loadings,
      method=method)
}
