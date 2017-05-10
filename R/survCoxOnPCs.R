survCoxOnPCs <- function(genes, expr, pcNum, annotations, formula,
                          pc2class=TRUE, robust=FALSE, shrink=FALSE, cliques=NULL) {

  expr <- expr[genes,, drop=FALSE]

  if (NROW(expr) == 0) {
    return(NULL)
  }

  expr <- t(expr)

  pcs <- computePCs(expr, pcNum , robust, shrink, cliques)

  apply(pcs, 2, survCoxOnPC, annotations=annotations, formula=formula, pc2class=pc2class)

}

survCoxOnPC <- function(pc, annotations, formula, pc2class) {

  if (pc2class){
    pc <- medianSurv(pc, 3)$class
  }

  coxObj <- data.frame(pc=pc, annotations)
  survivalcox(coxObj, formula)

}

survCoxOnAllPCs <- function(genes, expr, perc=0.8, annotations, pc2class=TRUE, robust=FALSE,
                            shrink=FALSE, cliques=NULL) {
  expr <- expr[genes,, drop=FALSE]

  if (NROW(expr) == 0) {
    return(NULL)
  }
  expr <- t(expr)
  pcs <- computePCs(expr, robust=robust, shrink=shrink, cliques=cliques)
  chosen <- choosePCS(pcs, perc)
  comps <- paste(chosen, collapse ="+")
  formula = as.formula(paste("Surv(days, status) ~", comps, sep=" "))
  pcsOnly <- pcs$x
  coxObj <- data.frame(pcs$x[,chosen,drop=F], annotations)
  if (length(chosen)==1){
    return(survivalcox(coxObj, formula))
  }
  fullsurvivalcox(coxObj, formula)
}
