\name{getMostCorrelatedGenes}
\alias{getMostCorrelatedGenes}

\title{Extract genes that are related to PCs associated with survival}
\description{Given a survCliques object the function extract those genes that are correlated the most with PCs that has been labeled as significant under a certain threshold.}

\usage{getMostCorrelatedGenes(scObj, thr=0.05, n=5)}

\arguments{
  \item{scObj}{an object survCliques}
  \item{thr}{thrshold to consider a clique significant. This threshold is used also for the significance of the zscores}
  \item{n}{retrun up to n top correlated genes}
}
\details{Function to reveal those genes that are more important in the survival process.}
\value{a list. All significant cliques. Correlations are expressed in gene heaviness of the loadings. Id exp is defined, also a correlation between significan PCs and gene profile is computed.}
\author{Paolo Martini}
