\name{distTable}
\alias{distTable}

\title{Distribution Table}

\description{
The function create a distribution table to be printed in Latex format.q
}

\usage{
distTable(df=NULL, df2=NULL, digits=4, dist=c("t", "chisq", "f", "normal"),
          size=0.05, CDF=FALSE, float=TRUE, label=NULL)
}

\arguments{

  \item{df}{A vector of degrees of freedom.}

  \item{df2}{A second vector of degrees of freedom for the f-distribution.}
  
  \item{digits}{How many digits to print?}

  \item{dist}{The distribution table to be printed.}

  \item{size}{The size for F-tests.}

  \item{CDF}{Should the normal table show CDF values? By default, it
  shows the probability of being between 0 and z.}

  \item{float}{Should the Latex table be a floating table?}

  \item{label}{The label of referencing the table in Latex. By default,
    it is \code{dist} followed by "Table".}
}
\value{
  It prints a latex table.
}

\examples{

\dontrun{distTable()}

}

