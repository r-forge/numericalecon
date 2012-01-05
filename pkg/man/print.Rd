\name{print}
\alias{print.Quadra}
\alias{print.zeros}
\alias{print.solveP.Quadra}
\title{Print methods for Quadra}
\description{
Print results of different methods applied to Quadra Objects
}
\usage{
\method{print}{Quadra}(x, ...)
\method{print}{zeros}(x, ...)
\method{print}{solveP.Quadra}(x, ...)
}
\arguments{
\item{x}{Onject of class \code{Quadra}.}
\item{...}{Argument for other type of objects}
}



\examples{
P1 <- Quadra(1,2,3)
print(P1)
print(zeros(P1))
print(solveP(P1))
}












