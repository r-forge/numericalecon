\name{inferenceReg}
\alias{inferenceReg}
\alias{testRegCoef}
\alias{ciRegCoef}
\alias{testLinRegCoefs}
\alias{testRegF}

\title{Object generators for inference solution in regression models}

\description{
The following are all solution generators for different
types of inference questions based on linear regression. They return an
object that can be converted into Latex though its print method.
}

\usage{

testRegCoef(object, h0, which, size=.05, alter=c("diff","greater","less"),
            digits=4, assume=c("Normal","nonNormal"),
            beta=NULL, se = NULL, n=NULL, ncoef=NULL, dround=NULL)

ciRegCoef(object, size=0.05,  which, digits=4,
          assume=c("Normal", "nonNormal"),
          beta=NULL, se=NULL, n=NULL, ncoef=NULL, dround=NULL)

testLinRegCoefs(object, R, q,  size=.05, alter=c("diff","greater","less"),
                digits=4, assume=c("Normal","nonNormal"),
                beta=NULL, vcoef=NULL, n=NULL, dround=NULL)

testRegF(object, objectr=NULL, size=.05, digits=4,
         assume = c("Normal", "nonNormal"), dround=NULL,
         Rsq=NULL, SSRu=NULL, SSRr=NULL, nrest=NULL, n=NULL, ncoef=NULL)
}

\arguments{

  \item{object,objectr}{Object of class \code{\link[stats]{lm}}. For
  F-tests, \code{objectr} is the \code{lm} object of the restricted
  model.}

  \item{h0}{The null hypothesis.}

  \item{which}{An integer specifying which coefficient is being tested.}
  
  \item{size}{The size of the test or one minus the confidence level
    when confidence intervals are computed.}

  \item{alter}{The alternative hypothesis.}

  \item{assume}{What is assumed about the distribution of the data. It
  affects how the critical values are obtained. The t or F distribution
  is used under normality, and the asymptotic N(0,1) or chi-square is
  used otherwise.}

  \item{digits}{How many digits to print?}

  \item{dround}{If we want to round the numbers before doing inference,
    \code{dround} is the number of decimals. This is particularly useful
    to generate solutions for exams in which coefficients, means and
    standard errors are rounded.}

  \item{beta,se,n,ncoef}{Alternatively, the coefficient estimate(s),
    the standard error(s), the sample size(s) and the number of coefficients
    in the regression model can be provided, instead of the \code{lm}
    object. The argument \code{se} can also be used to pass a vector
    of robust standard errors. By default, the non-robust one is used.}

  \item{R,q}{For test of linear combinations of coefficients, \code{R}
    is the hypothesis matrix with one row and \code{q} is the right-hand
    side.}

  \item{vcoef}{The alternative estimated covariance matrix. It is needed
    for robust testings}

  \item{Rsq}{The R-squared of the unrestricted model. If provided, an
  F-test for the significance of all coefficients but the intercept is
  performed.}
    
  \item{SSRu,SSRr,nrest}{The unrestricted SSR, the restricted SSR and
    the number of restrictions being tested. These values are needed if
    not \code{lm} objects are provided.}
}
  
\value{
  It returns an object of class \code{metricsSol}.
}

\examples{

## Test on one coefficient

data(wage1, package="wooldridge")
fit <- lm(wage~educ+exper, wage1)
testRegCoef(fit, h0=0.5, which=2, size=0.10, alter="greater")

## See vignette for more examples

}

