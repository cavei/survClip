\name{cliqueSurvivalTest}
\alias{cliqueSurvivalTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform survival test on all the cliques of the graph.}
\description{Given a pathway and a matrix with annotation pData compute survival test.}
\usage{cliqueSurvivalTest.Rd(expr, graph, survAnnot, pcNum=1, formula="Surv(days, status) ~ pc", pc2class=TRUE, root=NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{expression matrix.}
  \item{graph}{a graphNEL object.}
  \item{survAnnot}{survivalAnnot object.}
  \item{pcNum}{Numberof pc to consider. Default 1.}
  \item{formula}{the survival formula.}
  \item{pc2class}{if True confert the continous values in discrete.}
  \item{root}{The root for graph decomposition. If NULL gRbase gues it.}
}
\details{Survival test is made according to formula.}
\value{
  \item{alphas}{alphas for each cliques.}
  \item{cliques}{cliques definition}
}
\references{}
\author{Paolo Martini}
\seealso{\code{\link{survClip}}}
\examples{
  print("example")
}
