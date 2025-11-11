\name{inference}
\alias{inference}
\alias{testMean}
\alias{ciMean}
\alias{testDiffMeans}
\alias{ciDiffMeans}
\alias{testVar}
\alias{ciVar}
\alias{testDiffVar}

\title{Object generators for inference solution}

\description{
The following are all solution generators for different types of
inference questions. They return an object that can be converted into
Latex though its print method.
}

\usage{
testMean(x, h0, size=.05, alter=c("diff","greater","less"),
      assume=c("Normal","nonNormal"), digits=4, dround=NULL,
      xbar=NULL, se=NULL, n=NULL)

ciMean(x, size=0.05, assume=c("Normal","nonNormal"), digits=4, dround=NULL,
       xbar=NULL, se=NULL, n=NULL)

testDiffMeans(x1, x2, h0, size=.05, alter=c("diff","greater","less"),
              assume=c("Normal","nonNormal"), assumev=c("same", "diff"),
              digits=4, dround=NULL, xbar=NULL, se=NULL, n=NULL)

ciDiffMeans(x1, x2, size=.05, assume=c("Normal","nonNormal"),
            assumev=c("same", "diff"),
            digits=4, dround=NULL, xbar=NULL, se=NULL, n=NULL)

testVar(x, h0, size=.05, alter=c("diff","greater","less"),
        assume=c("Normal","nonNormal"), digits=4, dround=NULL,
        se=NULL, n=NULL)

ciVar(x, size=0.05, assume=c("Normal","nonNormal"), digits=4, dround=NULL,
      se=NULL, n=NULL)

testDiffVar(x1, x2, size=.05, 
            assume=c("Normal","nonNormal"), digits=4, dround=NULL,
            se=NULL, n=NULL)
}

\arguments{

  \item{x, x1, x2}{Numeric vectors}

  \item{h0}{The null hypothesis.}
  
  \item{size}{The size of the test or one minus the confidence level
    when confidence intervals are computed.}

  \item{alter}{The alternative hypothesis.}

  \item{assume}{What is assumed about the distribution of the data. It
  affects how the critical values are obtained. The t or F distribution
  is used under normality, and the asymptotic N(0,1) or chi-square is
  used otherwise.}

  \item{assumev}{Do we assume that the variance of the two samples are
  different or that they are same?}

  \item{digits}{How many digits to print?}

  \item{dround}{If we want to round the numbers before doing inference,
  \code{dround} is the number of decimals. This is particularly useful
  to generate solutions for exams in which coefficients, means and
  standard errors are rounded.}

  \item{xbar,se,n}{Alternatively, the sample mean(s), standard deviation(s)
    and sample size(s) can be provided, instead of numeric vector(s).
  }
}
\value{
  It returns an object of class \code{metricsSol}.
}

\examples{

## Test on the mean from the data

set.seed(112233)
x1 <- rnorm(30)
obj <- testMean(x1, 0)

## Test on the mean without data

obj <- testMean(h0=0, xbar=0.1, se=0.6, n=35)

## See vignette for more examples

}

