#' Perform survival test on all the cliques of the graph
#' 
#' This function performs survival test on given pathway using a matrix with survival annotation
#' 
#' Survival test is made according to survFormula. With "regular" method, a regular PCA analysis is used to compute PCs. With "sparse" method, a penalized regression is used for the estimation of PCs (as implemented in elasticnet)
#' 
#' @inheritParams pathwaySurvivalTest
#' 
#' @return A survCliques object
#' 
#' @seealso \code{\link{pathwaySurvivalTest}}, \code{\link{getTopLoadGenes}}
#' 
#' @examples
#' if (require(graphite)) {
#'   data(exp)
#'   data(survAnnot)
#'   data(graph)
#'   row.names(exp) <- paste0("ENTREZID:", row.names(exp))
#'   genes <- intersect(graph::nodes(graph), row.names(exp))
#'   graph <- graph::subGraph(genes, graph)
#'   expr <- exp[genes, , drop=FALSE]
#'   cliqueSurvivalTest(expr, survAnnot, graph, maxPCs=2)
#' }
#' 
#' @importFrom houseOfClipUtility extractCliquesFromDag
#' @importFrom graph nodes subGraph
#' @importFrom survival Surv
#' @importFrom methods new
#' 
#' @export
cliqueSurvivalTest <- function(expr, survAnnot, graph, pcsSurvCoxMethod=c("regular", "sparse"), alwaysShrink=FALSE, maxPCs=10, survFormula = "Surv(days, status) ~") {
  if (!is.data.frame(survAnnot)){
    stop("'annotations' must be a 'data.frame' object.")
  }
  
  pcsSurvCoxMethod <- pcsSurvCoxMethod[1]
  if (pcsSurvCoxMethod=="topological") {
    stop("topological method not supported for cliques.")
  }
  
  genes <- graph::nodes(graph)
  genes <- intersect(genes, row.names(expr))

  # decide if we want to stop or warn in no gene are found
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  graph <- graph::subGraph(genes, graph)
  expr <- expr[genes,, drop=FALSE]

  # clipper Function to import
  cliques <- houseOfClipUtility::extractCliquesFromDag(graph)
  results <- lapply(cliques, pcsSurvCox, expr=expr, annotations=survAnnot, method=pcsSurvCoxMethod, shrink=alwaysShrink, maxPCs=maxPCs, survFormula = survFormula)
  alphas  <- sapply(results, function(x) x$pvalue)
  zlist   <- lapply(results, function(x) x$zlist)
  cld     <- lapply(results, function(x) x$loadings)
  coxObjs <- lapply(results, function(x) x$coxObj)
  exprs   <- lapply(cliques, function(cls) {expr[cls, , drop=F]})
  
  names(alphas) <- NULL
  new("survCliques", alphas=alphas, zlist=zlist, cliques=cliques, coxObjs=coxObjs, cliquesLoadings=cld, cliquesExpr=exprs)
}
