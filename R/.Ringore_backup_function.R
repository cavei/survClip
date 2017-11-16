.DpathwaySurvivalTest <- function(expr, survAnnot, graph, pcsSurvCoxMethod=c("regular", "topological", "sparse"),
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
  
  set.seed(1234)
  gtpvalue          <- globaltest::p.value(globaltest::gt(Surv(days, events==1), alternative=t(expr)))
  set.seed(1234)
  pvalue            <- pcsSurvCox(genes, expr, survAnnot, method="regular", shrink=FALSE, cliques=NULL, maxPCs=maxPCs, survFormula = survFormula)
  set.seed(1234)
  pvalueShinkNoTopo <- pcsSurvCox(genes, expr, survAnnot, method="regular", shrink=TRUE, cliques=NULL, maxPCs=maxPCs, survFormula = survFormula)
  set.seed(1234)
  pvalueTopo        <- pcsSurvCox(genes, expr, survAnnot, method="topological", shrink=FALSE, cliques=cliques, maxPCs=maxPCs, survFormula = survFormula)
  set.seed(1234)
  pvalueTopoShrink  <- pcsSurvCox(genes, expr, survAnnot, method="topological", shrink=TRUE, cliques=cliques, maxPCs=maxPCs, survFormula = survFormula)
  
  new("survPath",
      pvalues = list(gtPvalue=gtpvalue, regPvalue=pvalue, regShrinkPvalue=pvalueShinkNoTopo, topoPvalue=pvalueTopo, topoShrinkPvalue=pvalueTopoShrink),
      method=method)
}

