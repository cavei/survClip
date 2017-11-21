\name{getTopLoadGenes}
\alias{getTopLoadGenes}

\title{Extract the relevant genes associated with survival}
\description{Given a survCliques object the function extracts those genes that are the most influent in the PCs identified as significant with a certain threshold.}

\usage{getTopLoadGenes(scObj, thr=0.05, n=5, loadThr=0.6)}

\arguments{
  \item{scObj}{an object survCliques}
  \item{thr}{threshold to consider a clique as significant. This threshold is used also for the significance of the zscores in zlist}
  \item{n}{return up to n top relevant genes}
  \item{loadThr}{filter loadings according to 'loadThr' absolute value}
}
\details{Function to reveal those genes that are more relevant in the
  survival process. The relevance of a gene is based on PC loadings.}

\value{a data.frame organized as follows: 
  \enumerate{
    \item{feature}{gene names}
    \item{clId}{clique id}
    \item{geneLoad}{gene loading}
    \item{whichPC}{the significant PC where the gene is relevant}
  }
  All significant cliques are represented. The importance of the genes is
  expressed by its loading.
}

\author{Paolo Martini}

\seealso{
\code{\link{cliqueSurvivalTest}}
}


\examples{
if (require(graphite)) {
  data(exp)
  data(survAnnot)
  data(graph)
  row.names(exp) <- paste0("ENTREZID:", row.names(exp))
  genes <- intersect(graph::nodes(graph), row.names(exp))
  graph <- graph::subGraph(genes, graph)
  expr <- exp[genes, , drop=FALSE]
  cliqueTest <- cliqueSurvivalTest(expr, survAnnot, graph, maxPCs=2)
  getTopLoadGenes(cliqueTest)
}
}