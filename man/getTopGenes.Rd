\name{getTopGenes}
\alias{getTopGenes}

\title{Extract genes that are related to PCs associated with survival}
\description{Given a survCliques object the function extract those genes that are correlated the most with PCs that has been labeled as significant under a certain threshold.}

\usage{getTopGenes(scObj, thr=0.05, n=5, corThr=0.6)}

\arguments{
  \item{scObj}{an object survCliques}
  \item{thr}{threshold to consider a clique significant. This threshold is used also for the significance of the zscores}
  \item{n}{retrun up to n top correlated genes}
  \item{corThr}{filter correlation/loadings according to this absolute value}
}
\details{Function to reveal those genes that are more important in the survival process.}
\value{a list. All significant cliques. Correlations are expressed in gene heaviness of the loadings. Id exp is defined, also a correlation between significan PCs and gene profile is computed.}
\author{Paolo Martini}
