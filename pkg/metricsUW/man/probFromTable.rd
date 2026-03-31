\name{probFromTable}
\alias{probFromTable}

\title{Probability from normal tables}

\description{
The function generates a Latex equation with the detailed steps to
calculate probability of areas using a printed normal table.
}

\usage{
probFromTable(a, b=NULL, type=c("less","greater"), digits=3,
              mu=0, sig2=1, varName="X", approx=FALSE, print=TRUE,
              CDF=FALSE)
}

\arguments{

  \item{a}{One of the boundary of the interval. See details.}

  \item{b}{One of the boundary of the interval. See details}

  \item{digits}{The number of digits shown in the normal table}

  \item{type}{When \code{b} is NULL, do we want the probability of being
  less or greater than \code{a}}
  
  \item{mu}{The expected value of the normal distribution for the
  unscaled variable.}

  \item{sig2}{The variance of the normal distribution for the
  unscaled variable.}

  \item{varName}{The name of the variable to be printed in the solution.}

  \item{approx}{Is the normal distribution an approximation?}

  \item{print}{Should the table be printed? If FALSE, the probabilities
    are returned.}

  \item{CDF}{Is the normal table showing the CDF? The alternative is a
  table that shows the probability of being between 0 and z.}
 
}

\details{
When \code{b} is NULL, it returns the probability that X is less than
\code{a} when \code{type="less"} and the probability of being greater
than \code{a} when \code{type="greater"}. When \code{b} is not NULL, it
returns the probability of being between \code{a} and \code{b}.
}

\value{
It returns an equation in Latex format or a numeric probability.
}

\examples{

\dontrun{probFromTable(-10, 8, mu=6, sig2=25, digits=3)}

}

