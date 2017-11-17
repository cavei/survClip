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

