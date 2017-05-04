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

computePCs <- function(exp, npc, robust, shrink=FALSE, cliques=NULL) {
  if (is.null(cliques)) {
    pcs <- hcomputePCs(exp, npc, robust)
  } else {
    pcs <- hcomputePCsTopo(exp, npc, shrink, cliques)
  }
  pcs
}
