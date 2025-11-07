\name{addQuadra}
\alias{addQuadra}
\alias{\%+\%}
\title{Sum of two polynomials}
\description{
It sums two objects of class Quadra
}
\usage{
addQuadra(Q1,Q2)
Q1\%+\%Q2
}
\arguments{
\item{Q1}{An object of class \code{Quadra}.}
\item{Q2}{An object of class \code{Quadra}.}
}

\examples{
P1 <- Quadra(1,2,3)
P2 <- Quadra(.1,-2,1.3)
P1%+%P2
}












