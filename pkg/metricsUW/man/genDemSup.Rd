\name{genDemSup}
\alias{genDemSup}
\alias{print.DemSup}
\alias{plot.DemSup}

\title{Demand and supply generator}

\description{
Functions to simulate data from a supply and demand system with one
exogenous variable in each function, plot the system and print 
the equations.
}

\usage{

genDemSup(n, a=c(20, 1, -2), d=c(100,-1,2,0), sige=16,
          sigeta=16, sigX=c(5,5), sigZ=5)

\method{print}{DemSup}(x, ...)

\method{plot}{DemSup}(x, y=NULL, nCurves=NULL, nPoints=NULL, fac=c(0.9,1.1), ...)
}

\arguments{

  \item{x}{An object of class \code{DemSup}.}

  \item{n}{The sample size.}

  \item{a}{The coefficient of the supply function. The first is the
  intercept, the second is the slope of the supply and the third is the
  coefficient of the exogenous variable Z.}

  \item{d}{The coefficient of the demand function. The first is the
  intercept, the second is the slope of the demand and the third is the
  coefficient of the exogenous variable X.}

  \item{sige,sigeta,sigX,sigZ}{The variance of the error in the supply,
  the error in the demand, the exogenous variable X and the exogenous
  variable Z.}

  \item{nCurves}{How many supply and demand curves to plot. By default,
    all curves are plotted.}

  \item{nPoints}{How many data points to plot. By default, they are all
    plotted.}

  \item{fac}{It determine the range for each demand and supply
  function. The function goes from \code{fac[1]} times the observed
  price to \code{fac[2]} times the observed price.}

  \item{y}{Not used at the moment.}

  \item{...}{Arguments to pass to other methods. Currently not used.}
}

\value{
\code{genDemSup} returns an object of class \code{DemSup}. The object is
a list that contains the following:

\item{dat}{A \code{data.frame} with the quantity Q, the price P, the
  exogenous variable from the demand X, the exogenous variable from the
  supply Z, the error from the demand eta and the error from the supply
  e.}

\item{a}{The parameters of the supply.}

\item{d}{The parameters of the demand.}

\item{sig}{The variance of the error terms.}

\item{sigExo}{The variance of the exogenous variables.}

}

\examples{

set.seed(112233)
sys <- genDemSup(200)
print(sys)
plot(sys, nCurves=10)

}

