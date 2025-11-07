\name{plot}
\alias{plot.Quadra}
\title{Plot second order polynomial}
\description{
It prints the polynomial and shows the min or max and the zeros.
}
\usage{
\method{plot}{Quadra}(x,from=NULL,to=NULL, ...)
}

\arguments{
  \item{x}{Object of class Quadra.}
  \item{from}{min(x) on the plot}
  \item{to}{max(x) on the plot}
\item{...}{Argument for other type of objects}
  }

\examples{

P1 <- Quadra(1,2,3)
plot(P1)
}
