cliqueSurvivalTest <- function(expr, survAnnot, graph, pcsSurvCoxMethod=c("regular", "sparse"), alwaysShrink=FALSE, maxPCs=10, survFormula = "Surv(days, status) ~") {
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
  cliques <- clipper:::extractCliquesFromDag(graph)
  results <- lapply(cliques, pcsSurvCox, expr=expr, annotations=survAnnot, method=pcsSurvCoxMethod, shrink=alwaysShrink, maxPCs=maxPCs, survFormula = survFormula)
  alphas  <- sapply(results, function(x) x$pvalue)
  zlist   <- lapply(results, function(x) x$zlist)
  cld     <- lapply(results, function(x) x$loadings)
  coxObjs <- lapply(results, function(x) x$coxObj)
  exprs   <- lapply(cliques, function(cls) {expr[cls, , drop=F]})
  
  names(alphas) <- NULL
  new("survCliques", alphas=alphas, zlist=zlist, cliques=cliques, coxObjs=coxObjs, cliquesLoadings=cld, cliquesExpr=exprs)
}
