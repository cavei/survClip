\name{getTopLoadGenes}
\alias{getTopLoadGenes}

\title{Extract the genes that are related to the PCs founded associated with survival}
\description{Given a survCliques object the function extracts those genes that are the most influent in the PCs that has been labeled as significant under a certain threshold.}

\usage{getTopLoadGenes(scObj, thr=0.05, n=5, loadThr=0.6)}

\arguments{
  \item{scObj}{an object survCliques}
  \item{thr}{threshold to consider a clique significant. This threshold is used also for the significance of the zscores in zlist}
  \item{n}{return up to n top correlated genes}
  \item{loadThr}{filter loadings according to this absolute value}
}
\details{Function to reveal those genes that are more important in the survival process. The importance of a gene is its weight in the PC loading.}
\value{a list. All significant cliques. The importance of the genes is expressed by its load.}
\author{Paolo Martini}
\examples{
if (require(graphite)) {
  data(exp)
  data(survAnnot)
  row.names(exp) <- paste0("ENTREZID:", row.names(exp))
  k <- pathways("hsapiens", "kegg")
  p <- convertIdentifiers(k[["Pathways in cancer"]], "entrez")
  graph <- pathwayGraph(p)
  genes <- intersect(graph::nodes(graph), row.names(exp))
  graph <- graph::subGraph(genes, graph)
  expr <- exp[genes, , drop=FALSE]
  cliqueTest <- cliqueSurvivalTest(expr, survAnnot, graph, maxPCs=2)
  getTopLoadGenes(cliqueTest)
}
}