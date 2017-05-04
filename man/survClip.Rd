\name{survClip}
\alias{survClip}
\title{Dissect the pathway to find the path with the greatest
  association with patients survival.}
\description{Basing on survival cliques analysis, this
  function identifies the paths that are mostly related with prognosis.}
\usage{survClip(expr, survAnnot, graph, pcNum=1, formula="Surv(days, status) ~ pc",
  pc2class=TRUE, nperm=100, roots=NULL, trZero=0.001, signThr=0.05, maxGap=1, dropNULL=FALSE)}
\arguments{
  \item{expr}{an expression  matrix or ExpressionSet with colnames for
    samples and row name for genes.}
  \item{survAnnot}{A data.frame with days (mandatory) and other classes for survival.}
  \item{graph}{a \code{graphNEL} object.}
  \item{pcNum}{the number of PCs to consider. Default 1.}
  \item{formula}{the Surv formula to be applied.}
  \item{pc2class}{If true convert continuos malues in classes. Default TRUE.}
  \item{nperm}{number of permutations. Default = 100.}
  \item{roots}{vector of possible root genes. Default NULL means the choice is made by gRbase.}
  \item{trZero}{lowest pvalue detectable. This threshold avoids that
    -log(p) goes infinite.}
  \item{signThr}{significance threshold for clique pvalues.}
  \item{maxGap}{allow up to maxGap gaps in the best path
    computation. Default = 1.}
  \item{dropNULL}{Default FALSE. Remove NULL root from the results list.}
}

\details{
  The method applied is a cox regression following formula.
}

\value{A matrix with a row for each paths. Columns are organized
  as follows:
  \enumerate{
    \item Index of the starting clique
    \item Index of the ending clique
    \item Index of the clique where the maximum value is reached
    \item Length of the path
    \item Maximum score of the path
    \item Average score along the path
    \item Percentage of path activation
    \item Impact of the path on the entire pathway
    \item Cliques involved and significant
    \item Cliques forming the path
    \item Genes forming the significant cliques
    \item Genes forming the path
  }
}

