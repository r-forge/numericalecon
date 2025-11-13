\name{printReg}
\alias{printReg}
\alias{printReg.lm}
\alias{printReg.tsls}
\alias{printReg.formula}

\title{Printing Estimated Models}

\description{
The function generates a Latex equation from a model estimated by
\code{\link[stats]{lm}} or \code{\link[stats]{glm}}.
}

\usage{
\method{printReg}{lm}(obj, digits=4, format="f", maxpl=5,
    adjrsq=FALSE, se=NULL, stars=FALSE, dist=c("t","n"), label=NULL,
    omit=NULL, ...)

\method{printReg}{tsls}(obj,digits=4, format='f', maxpl=5, 
    se=NULL, stars=FALSE, label=NULL, omit=NULL,
    strength=FALSE, Jtest=FALSE, Exo=FALSE, ...)

\method{printReg}{formula}(obj, maxpl=5, label=NULL, ...)
}

\arguments{

  \item{obj}{An object used to print the regression equation.}

  \item{digits}{The number of digits to be printed (See
          \code{\link{formatC}}).}
  
  \item{format}{Printing format for coefficient estimates and their
          standard errors.(See \code{\link{formatC}})}

  \item{maxpl}{The maximum number of variables to be printed per line.}

  \item{adjrsq}{Should we print the adjusted R-square?}

  \item{se}{Alternative vector of standard errors. By default, the one
  from \code{\link[stats]{vcov}} is used.}
  
  \item{stars}{Should we add stars next to significant coefficients?}
  
  \item{dist}{Distribution used for p-values. By default, the
  t-distribution is used. Alternatively, we can replace it by the
  asymptotic N(0,1).}

  \item{label}{Label for latex equations in case a reference is
  needed. By default, it is null and the equation is not numbered.}

  \item{omit}{This argument is to specify variables to omit from the
    printed equation. It is a character vector, and the function omits
    any variables that contain the characters detected by
    \code{\link{grep}}. When specified, the equation ends by specifying
    how many variables are omitted.}

  \item{strength, Jtest, Exo}{Should the first stage F-test, J-test and
  the list of exogenous variables be printed?}

  \item{...}{Other argument to pass to \code{\link{format}}}
}

\value{
It returns an equation in Latex format.
}

\examples{

## Estimated linear model

set.seed(112233)
x1 <- rnorm(30)
x2 <- rnorm(30)
y <- x1+x2+rnorm(30)
fit <- lm(y~x1+x2)
printReg(fit)

## formula to equation

printReg(y~x1*x2+z+I(z^2))
}

