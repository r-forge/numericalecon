\name{solveP}
\alias{solveP.Quadra}
\alias{solveP}
\title{Roots of a polynomial}
\description{
Computes the roots of a second order polynomial
}
\usage{
\method{solveP}{Quadra}(obj, ...)
}
\arguments{
\item{obj}{An object of class \code{Quadra}.}
\item{...}{Other arguments when \code{solve} is applied to another class object}
}

\examples{
P1 <- Quadra(1,2,3)
solveP(P1)
}












