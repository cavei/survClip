hcomputePCs <- function(exp, npc, robust){
  if (NCOL(exp) < npc)
    npc <- NCOL(exp)

  if(robust){
    rrcov::PcaHubert(exp,k=npc)@scores
  } else {
    pc1 <- stats::prcomp(exp)$x[,1:npc]
    as.matrix(pc1,ncol=npc)
  }
}

hcomputePCsTopo <- function(exp, npc, shrink, cliques){
  if (NCOL(exp) < npc)
    npc <- NCOL(exp)

  covmat <- clipper:::estimateExprCov(exp, shrink)
  cliquesIdx <- lapply(cliques, function(c) match(c, row.names(covmat)))
  scovmat <- qpgraph::qpIPF(covmat, cliquesIdx)

  pcCov<-base::eigen(scovmat)
  scalee <- scale(exp)
  eigenvector <- pcCov$vectors
  scores <- scalee%*%eigenvector
  scores[, seq_len(npc), drop=FALSE]
}

compPCs <- function(exp, robust){
  if(robust){
    pca <- rrcov::PcaHubert(exp)
    sdev <- apply(pca@scores,2,sd)
    return(list(x=pca@scores, sdev=sdev))
  }
  pca <- stats::prcomp(exp)
  return(list(x=pca$x, sdev=pca$sdev))
}

compPCsTopo <- function(exp, shrink, cliques){
  covmat <- clipper:::estimateExprCov(exp, shrink)
  covmat <- makePositiveDefinite(covmat)$m1
  cliquesIdx <- lapply(cliques, function(c) match(c, row.names(covmat)))
  scovmat <- qpgraph::qpIPF(covmat, cliquesIdx)
  pcCov<-base::eigen(scovmat)
  scalee <- scale(exp)
  eigenvector <- pcCov$vectors
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector[, seq_len(npc), drop=F]
  colnames(scores) <- paste0("PC", seq_len(npc))
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd))
}

choosePCS<-function(pcs, variability) {
  percentVar <- pcs$sdev^2 / sum(pcs$sdev^2)
  names(percentVar) <- colnames(pcs$x)
  selection <- !cumsum(percentVar)>=variability
  if (sum(selection) == 0) {
    return("PC1")
  }
  return(names(which(selection)))
}

# computePCs <- function(exp, npc, robust, shrink=FALSE, cliques=NULL) {
#   if (is.null(cliques)) {
#     pcs <- hcomputePCs(exp, npc, robust)
#   } else {
#     pcs <- hcomputePCsTopo(exp, npc, shrink, cliques)
#   }
#   pcs
# }

computePCs <- function(exp, robust=FALSE, shrink=FALSE, cliques=NULL) {
  if (is.null(cliques))
    return(compPCs(exp, robust))
  pcs <- compPCsTopo(exp, shrink, cliques)
}
