chooseBestPC <- function(alpha){
  idx <- which.min(alpha)
  p <- min(alpha)
  return(c(p, idx))
}

cliqueSurvivalTest <- function(expr, graph, survAnnot, pcNum=1, formula="Surv(days, status) ~ pc", pc2class=TRUE, root=NULL) {
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
  alphas <- lapply(cliques, survCoxOnPCs, expr=expr, pcNum=pcNum, annotations=survAnnot, formula=formula, pc2class=pc2class)
  alphas <- t(sapply(alphas, chooseBestPC))

  list(alpha=alphas[,1], pc=alphas[,2], cliques=cliques)
}
