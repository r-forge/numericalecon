\name{zeros}
\alias{zeros.Quadra}
\alias{zeros}
\title{Zeros of a second order polynomial}
\description{
Computes the zeros of polynomials of class Quadra
}
\usage{
\method{zeros}{Quadra}(object, ...)
}
\arguments{
\item{object}{An object os class \code{Quadra}}
\item{...}{Other arguments when \code{zeros} is applied to an other class object}
}

\examples{

P1 <- Quadra(1,2,3)
zeros(P1)

}

