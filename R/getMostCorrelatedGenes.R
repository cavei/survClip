correlateGeneToPC <- function(pcs, loadings, n) {
  lapply(pcs, function(pc) {
    ld.pc <- loadings[, pc, drop=F]
    fout <- row.names(ld.pc)[which(ld.pc == 0)]
    genes <- row.names(ld.pc)[head(order(abs(ld.pc), decreasing = TRUE), n)]
    ld.pc[setdiff(genes, fout), , drop=F]
  })
}

corPCandGenes <- function(pcs, coxObj, geneExp, n) {
  lapply(pcs, function(pc) {
    pc.value <- coxObj[, pc, drop=F]
    correlations <- t(cor(pc.value, t(geneExp)))
    selected <- head(order(abs(correlations), decreasing = TRUE), n)
    correlations[selected, , drop=F]
  })
}

getMostCorrelatedGenes <- function(scObj, thr=0.05, n=5) {
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
    ldCor <- correlateGeneToPC(pcs, loadings, n)
    pcCor <- corPCandGenes(pcs, coxObj, exprs[[clId]], n)
    list(cliqueId=clId, loadingsBased=ldCor, correlationBased=pcCor)
  })
}
