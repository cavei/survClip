\name{makePositiveDefinite}
\alias{makePositiveDefinite}
\title{make positive and defined a matrix}
\description{for internal use only}
\usage{
  makePositiveDefinite(m1, m2 = NULL, m3 = NULL, threshold = 0.1)
}
\arguments{
  \item{m1}{matrix}
  \item{m2}{optional. Second matrix.}
  \item{m3}{optional. Third matrix.}
  \item{threshold}{threshold}
}
\details{for internal use only.}
\value{ a list with 
\item{m1}{matrix}
\item{m2}{matrix}
\item{m3}{matrix}
\item{correction}{logical. If a correction was made}
\item{value}{correction magnetude}
}
\author{
Paolo Martini
}