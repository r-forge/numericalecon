\name{print}
\alias{print.metricsSol}

\title{Print method for objects of class \code{metricsSol}}

\description{
It prints the solution from the object using Latex format. It is mainly
used inside markdown documents to generate nicely formatted solution.
}

\usage{
\method{print}{metricsSol}(x, addMess=NULL, ...)
}

\arguments{

  \item{x}{An object of class \code{metricsSol}.}

  \item{addMess}{Message to print at the end of the solution. I use it
    to add the mark distribution when I generate solutions for exams.}

  \item{...}{Arguments to pass to other methods. Currently not used.}
}

\value{
It print the type of inference and the detailed solution in Latex format.
}

\examples{

set.seed(112233)
x1 <- rnorm(30)
obj <- testMean(x1, 0)
print(obj, addMess="Not point if no justification")
}

