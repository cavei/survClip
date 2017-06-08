cliqueSurvivalTest <- function(expr, graph, survAnnot, root=NULL, pcsSurvCoxMethod=c("regular", "sparse"), alwaysShrink=FALSE, maxPCs=10) {
  if (!is.data.frame(survAnnot)){
    stop("'annotations' must be a 'data.frame' object.")
  }
  pcsSurvCoxMethod <- pcsSurvCoxMethod[1]
  if (pcsSurvCoxMethod=="topological") {
    stop("topological method not supported for cliques.")
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
  results <- lapply(cliques, pcsSurvCox, expr=expr, annotations=survAnnot, method=pcsSurvCoxMethod, shrink=alwaysShrink, maxPCs=maxPCs)
  alphas <- sapply(results, function(x) x$pvalue)
  zlist  <- lapply(results, function(x) x$zlist)
  names(alphas) <- NULL
  new("survCliques", alphas=alphas, zlist=zlist, cliques=cliques)
}
