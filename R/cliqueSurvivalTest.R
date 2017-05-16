chooseBestPC <- function(alpha){
  idx <- which.min(alpha)
  p <- min(alpha)
  return(c(p, idx))
}

cliqueSurvivalTest <- function(expr, graph, survAnnot, pcNum=1, perc=0.6, formula="Surv(days, status) ~ pc", pc2class=TRUE, root=NULL) {
  if (!is.data.frame(survAnnot)){
    stop("'annotations' must be a 'data.frame' object.")
  }

  genes <- nodes(graph)
  genes <- intersect(genes, row.names(expr))

  # decide if we want to stop or warn in no gene are found
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  graph <- graph::subGraph(genes, graph)
  expr <- expr[genes,, drop=FALSE]

  # clipper Function to import
  cliques <- clipper:::extractCliquesFromDag(graph, root=root)
  
  # survCoxOnAllPCs(genes, expr, perc=perc, annotations, pc2class=TRUE, robust=FALSE, shrink=FALSE, cliques=NULL)
  results <- lapply(cliques, survCoxOnAllPCs, expr=expr, perc=perc, annotations=survAnnot, pc2class=pc2class, robust=robust, shrink=FALSE, cliques=NULL)
  alphas <- sapply(results, function(x) x$pvalue)
  names(alphas) <- NULL
  list(alpha=alphas, full=results, cliques=cliques)
}
