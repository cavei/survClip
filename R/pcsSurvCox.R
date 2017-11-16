pcsSurvCox <- function(genes, expr, annotations, method=c("regular", "topological", "sparse"), shrink=FALSE,cliques=NULL, maxPCs=10,survFormula = "Surv(days, status) ~") {
  expr <- expr[genes,, drop=FALSE]
  
  if (NROW(expr) == 0) {
    return(NULL)
  }
  expr <- t(expr) ## check this

  if (NCOL(expr)!=1) {
    pcs <- computePCs(expr, shrink=shrink, method=method, cliques=cliques, maxPCs=maxPCs)
  } else {
    colnames(expr) <- "PC1"
    pcs <- list(x=expr, sdev=sd(expr), loadings=1)
  }
  
  comps <- paste(colnames(pcs$x), collapse ="+")
  formula = as.formula(paste(survFormula, comps, sep=" "))
  coxObj <- data.frame(pcs$x, annotations)
  scox <- survivalcox(coxObj, formula)
  scox$loadings <- pcs$loadings
  scox
}
