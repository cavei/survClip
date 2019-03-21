#' Survival analysis on the whole pathway
#' 
#' This function performs survival analysis on pathways. The analysis can be either topological or classical. The analysis is based on data reduction based on Principal Component Analysis and a Cox proportional hazard model on the most influent PCs
#' 
#' With "regular" method, a regular PCA analysis is used to compute PCs. With "topological" method, the covariance matrix is estimated using the topology of the pathway with IPS algorithm. With "sparse" method, a penalized regression is used for the estimation of PCs (as implemented in elasticnet). The max number of PCs used by the model is estimated by "estim_ncp" in FactoMineR. A maximum number of PCs can be fixed by the user. The minimum ot the two is chosen.
#' 
#' @param expr expression matrix
#' @param survAnnot a data frame for survival annotations specified according to the survFormula. The data frame must contain days and status
#' @param graph a graphNEL object for a graph
#' @param pcsSurvCoxMethod a method to perform PCA. Can be "regular", "topological", "sparse" for regular PCA, topological based PCA and sparse PCA, respectively. The latter one (sparse) is particularly suited for cliques only
#' @param alwaysShrink if TRUE, always shrink the covariance matrix. Deafult=FALSE
#' @param maxPCs maximum number of PCs used in the cox formula "Surv(days, status) ~ PC1.."
#' @param survFormula the formula used in Coxph analysis. Defaut="Surv(days, status) ~". Please note that the formula end with '~' meaning that PCs will be added
#' @param robust should be used the robust mode for cox
#' 
#' @return A survPath object
#' 
#' @seealso \code{\link{cliqueSurvivalTest}}
#' @examples
#' if (require(graphite)) {
#'   data(exp)
#'   data(survAnnot)
#'   data(graph)
#'   row.names(exp) <- paste0("ENTREZID:", row.names(exp))
#'   genes <- intersect(graph::nodes(graph), row.names(exp))
#'   graph <- graph::subGraph(genes, graph)
#'   expr <- exp[genes, , drop=FALSE]
#'   pathwaySurvivalTest(expr, survAnnot, graph, maxPCs=2)
#' }
#' 
#' @importFrom houseOfClipUtility extractCliquesFromDag
#' @importFrom graph nodes subGraph
#' @importFrom survival Surv
#' @importFrom methods new
#' 
#' @export
#' 
pathwaySurvivalTest <- function(expr, survAnnot, graph, pcsSurvCoxMethod=c("regular", "topological", "sparse"),
                                alwaysShrink=FALSE, maxPCs=10, survFormula = "Surv(days, status) ~",
                                robust=FALSE){
  genes <- graph::nodes(graph)
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

  cliques <- houseOfClipUtility::extractCliquesFromDag(graph)

  maxcliques <- max(sapply(cliques, length))
  shrink <- length(samples) < maxcliques | alwaysShrink

  days   <- survAnnot$days
  events <- survAnnot$status
  
  method = pcsSurvCoxMethod[1]
  
  res <- pcsSurvCox(genes, expr, survAnnot, method=method, shrink=alwaysShrink, cliques=cliques, maxPCs=maxPCs, survFormula = survFormula, robust=robust)
  new("survPath",
      pvalue = res$pvalue, zlist = res$zlist, coxObj = res$coxObj, loadings = res$loadings,
      method=method)
}
