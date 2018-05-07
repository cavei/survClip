#' Cox Analysis on Principal Components
#' 
#' Survival Analysis using Cox Model on Principal Components. For internal use only.
#' 
#' For internal use only.
#' 
#' @param genes vector of genes
#' @param expr whole expression matrix
#' @param annotations survival annotations
#' @param method PCA methods
#' @param shrink shrink covariance matrix
#' @param cliques when method=topological use this cliques to define topology in IPF
#' @param maxPCs maximum number of PC to consider
#' @param survFormula survival forumula to use
#' 
#' @return cox model object
#' 
#' @importFrom houseOfClipUtility computePCs
#' @importFrom stats as.formula sd na.omit
#' @importFrom survival coxph Surv 
#' 
#' @export
#' 
pcsSurvCox <- function(genes, expr, annotations, method=c("regular", "topological", "sparse"), shrink=FALSE,cliques=NULL, maxPCs=10,survFormula = "Surv(days, status) ~") {
  expr <- expr[genes,, drop=FALSE]
  
  if (NROW(expr) == 0) {
    return(NULL)
  }
  expr <- t(expr) ## check this

  if (NCOL(expr)!=1) {
    pcs <- houseOfClipUtility::computePCs(expr, shrink=shrink, method=method, cliques=cliques, maxPCs=maxPCs)
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
