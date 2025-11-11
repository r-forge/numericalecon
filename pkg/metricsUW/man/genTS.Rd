\name{genTS}
\alias{genTS}

\title{Functions to generate datasets}

\description{
The function generate datasets with known properties. Students can try
different techniques to identify the processes.
}

\usage{

genTS(n=200, frequency=12, start=c(1960,1), seed)

}

\arguments{

  \item{n}{The sample size}

  \item{frequency,start}{The frequency and starting data. See \code{\link[stats]{ts}}.}
  
  \item{seed}{An optional seed for the random generators}

}
  
\value{
  It returns an \code{mts} object with two time series, \code{X} and
  \code{Y}. One of them is a unit root process and the other is a trend
  stationary series. The true processes are saved in the attribute "which".
}

\examples{

dat <- genTS()
}

