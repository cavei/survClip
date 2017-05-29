checkIn <- function(expr, graph, root) {
  if (NROW(expr)==0){
    warning("Expression matrix has 0 rows.")
    return(NULL)
  }
  if (!identical(row.names(expr), graph::nodes(graph))){
    stop("Mismatch between graph nodes and expression row nams.")
  }
  if (!is.null(root)){
    if (!(root %in% row.names(expr)))
      stop(paste(root, ' is not a valid root', sep=""))
  }
  return(TRUE)
}

cleanClipperResults <- function(clipped) {
  clpprNames <- c("startIdx", "endIdx", "maxIdx", "lenght", "maxScore", "aScore", "activation", "impact",
    "involvedCliques", "cliqueOnPath", "involvedGenes", "pathGenes")

  if (!is.matrix(clipped))
    clipped <- as.matrix(clipped)

  clipped <- t(clipped)
  colnames(clipped) <- clpprNames
  clipped <- clipper:::addNames(clipped)
  clipped <- clipper:::removeNArows(clipped)
  if (NROW(clipped)==0)
    return(NULL)

  clipped <- as.data.frame(clipped, stringsAsFactors=FALSE)
  as.data.frame(clipped, stringsAsFactors=FALSE)
}

singleSurvivalClip <- function(root, expr, survAnnot, graph, pcNum, perc, formula, pc2class, nperm, trZero, signThr, maxGap, robust) {
  checkIn(expr, graph, root)
  ct <- cliqueSurvivalTest(expr, graph, survAnnot, pcNum, perc, formula=formula, pc2class, robust, root)
  if (is.null(ct)){
    return(NULL)
  }

  jtp <- clipper::getJunctionTreePaths(graph, root)
  if (is.null(jtp))
    return(NULL)

  clipped <- clipper:::runCoreClipper(ct, jtp, trZero, signThr, maxGap)
  cleanClipperResults(clipped)
}

chooseRoot <- function(allTests) {
  maxScore <- 0
  maxRoot <- NULL
  for (r in names(allTests)) {
    if (NROW(allTests[[r]])!=0){
      maxRelativeScore <- max(as.numeric(allTests[[r]][,5]))
      if (maxScore < maxRelativeScore) {
        maxScore <- maxRelativeScore
        maxRoot <- r
      }
    }
  }
  if (is.null(maxRoot))
    return(NULL)

  row.names(allTests[[maxRoot]]) <- paste(maxRoot, row.names(allTests[[maxRoot]]), sep="-")
  allTests[[maxRoot]]
}

survClip <- function(expr, survAnnot, graph, pcNum=1, perc=0.6, formula="Surv(days, status) ~ pc", pc2class=TRUE,
                     nperm=100, roots=NULL, trZero=0.001, signThr=0.05, maxGap=1, dropNULL=FALSE, robust=FALSE, shrinkForCLiques=FALSE){

  # pcNum=1; perc=0.6; formula="Surv(days; status) ~ pc"; pc2class=TRUE; nperm=100; roots=NULL; trZero=0.001; signThr=0.05; maxGap=1; dropNULL=FALSE
  
  if (is.null(roots) & dropNULL) {
    stop('You can not drop NULL as root when your only root is null.')
  }
  # Need a function for the check?
  if (NROW(expr)==0){
    warning("Expression matrix has 0 rows.")
    return(NULL)
  }

  expGenes <- row.names(expr)
  genes <- graph::nodes(graph)
  genes <- intersect(genes, expGenes)

  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  graph <- graph::subGraph(genes, graph)

  if (is.null(roots)) {
    rootNULL <- singleSurvivalClip(root=NULL, expr=expr, survAnnot=survAnnot, graph=graph, pcNum=pcNum, perc=perc, formula=formula,
                                   pc2class=pc2class, nperm=nperm, trZero=trZero, signThr=signThr, maxGap=maxGap, robust=robust, shrinkForCliques=shrinkForCliques)
    return(rootNULL)
  }

  allTests <- lapply(roots, singleSurvivalClip, expr=expr, survAnnot=survAnnot, graph=graph, pcNum=pcNum, perc=perc, formula=formula,
                    pc2class=pc2class, nperm=nperm, trZero=trZero, signThr=signThr, maxGap=maxGap, robust=robust, shrinkForCliques=shrinkForCliques)

  names(allTests) <- roots

  if (!dropNULL) {
    rootNULL <- singleSurvivalClip(root=NULL, expr=expr, survAnnot=survAnnot, graph=graph, pcNum=pcNum, perc=perc, formula=formula,
                                   pc2class=pc2class, nperm=nperm, trZero=trZero, signThr=signThr, maxGap=maxGap, robust=robust, shrinkForCliques=shrinkForCliques)
    rnull <- length(allTests) + 1
    allTests[[rnull]] <- rootNULL
    names(allTests) <- c(roots,"null")
  }

  chooseRoot(allTests)
}
