#' Extract the relevant genes associated with survival
#' 
#' Given a survCliques object the function extracts those genes that are the most influent in the PCs identified as significant with a certain threshold
#' 
#' Function to reveal those genes that are more relevant in the survival process. The relevance of a gene is based on PC loadings
#' 
#' @param scObj an object survCliques
#' @param thr threshold to consider a clique as significant. This threshold is used also for the significance of the zscores in zlist
#' @param n return up to n top relevant genes
#' @param loadThr filter loadings according to 'loadThr' absolute value
#' 
#' @return a data.frame organized as follows: 
#' \enumerate{
#'   \item{feature}{gene names}
#'   \item{clId}{clique id}
#'   \item{geneLoad}{gene loading}
#'   \item{whichPC}{the significant PC where the gene is relevant}
#' }
#' All significant cliques are represented. The importance of the genes is expressed by its loading.
#' 
#' @seealso \code{\link{cliqueSurvivalTest}}
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
#'   cliqueTest <- cliqueSurvivalTest(expr, survAnnot, graph, maxPCs=2)
#'   getTopLoadGenes(cliqueTest)
#' }
#' 
#' @importFrom checkmate assertClass
#' @importFrom utils head
#' 
#' @export
#' 
getTopLoadGenes <- function(scObj, thr=0.05, n=5, loadThr=0.6) {
  assertClass(scObj, "survCliques")
  idx <- which(scObj@alphas <= thr)
  if (length(idx) == 0)
    return(NULL)
  
  ld <- scObj@cliquesLoadings
  z  <- scObj@zlist
  coxObjs <- scObj@coxObjs
  exprs <- scObj@cliquesExpr
  
  corGenes <- lapply(idx, function(clId){
    loadings <- ld[[clId]]
    coxObj <- coxObjs[[clId]]
    pcs <- names(which(z[[clId]] <= thr))
    ldCors <- correlateGeneToPC(pcs, loadings, n, loadThr)
    if (length(ldCors)==0)
      return(c(clId, "NULL", "NULL"))
    
    form <- lapply(ldCors, function(ldCor) {
      rm <- cbind(clId, ldCor, pc=colnames(ldCor))
      row.names(rm) <- row.names(ldCor)
      colnames(rm)[2] <- "ld"
      rm
    })
    do.call(rbind, form)
  })
  mat <- do.call(rbind, corGenes)
  data.frame(feature=row.names(mat), clId=mat[,1], geneLoad=mat[,2], whichPC=mat[,3])
  
}

correlateGeneToPC <- function(pcs, loadings, n, thr) {
  lapply(pcs, function(pc) {
    ld.pc <- loadings[, pc, drop=F]
    fout <- row.names(ld.pc)[which(ld.pc == 0)]
    genes <- row.names(ld.pc)[head(order(abs(ld.pc), decreasing = TRUE), n)]
    selectLoad <- ld.pc[setdiff(genes, fout), , drop=F]
    selection <- selectLoad[,1] <= -1*thr | selectLoad >= thr
    selectLoad[selection, , drop=F]
  })
}

#' @importFrom stats cor
corPCandGenes <- function(pcs, coxObj, geneExp, n, thr) {
  lapply(pcs, function(pc) {
    pc.value <- coxObj[, pc, drop=F]
    correlations <- t(cor(pc.value, t(geneExp)))
    selected <- head(order(abs(correlations), decreasing = TRUE), n)
    fullC <- correlations[selected, , drop=F]
    selection <- fullC[,1] <= -1*thr | fullC[,1] >= thr
    fullC[selection, , drop=F]
  })
}

getTopGenes <- function(scObj, thr=0.05, n=5, corThr=0.6) {
  idx <- which(scObj@alphas <= thr)
  if (length(idx) == 0)
    return(NULL)
  
  ld <- scObj@cliquesLoadings
  z  <- scObj@zlist
  coxObjs <- scObj@coxObjs
  exprs <- scObj@cliquesExpr
  
  lapply(idx, function(clId){
    loadings <- ld[[clId]]
    coxObj <- coxObjs[[clId]]
    pcs <- names(which(z[[clId]] <= thr))
    ldCor <- correlateGeneToPC(pcs, loadings, n, corThr)
    pcCor <- corPCandGenes(pcs, coxObj, exprs[[clId]], n, corThr)
    list(cliqueId=clId, loadingsBased=ldCor, correlationBased=pcCor)
  })
}

