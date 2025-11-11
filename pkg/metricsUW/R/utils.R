## TS generator

genTS <- function(n=200, frequency=12, start=c(1960,1), seed)
{
    if (!missing(seed))
        set.seed(seed)    
    u1 <- arima.sim(n, model=list(ar=runif(1,.3,.95)))*5
    beta <- c(runif(1, 5, 15), sample(c(0, runif(1, -5,5)), 1))
    u2 <- beta[1] + beta[2]*u1 + arima.sim(n, model=list(ar=runif(1,.3,.95)))*5
    d1 <- c(runif(1, 0, 5),
            runif(1, -1,))
    tmp <- sample(c(0, ifelse(d1[2]<0, runif(1, -2*n*d1[2], -2*n*d1[2]*3),
                              runif(1, -2*d1[2], 0))),1)
    d1 <- c(d1[1], tmp, d1[2])
    T1 <- d1[1]+d1[2]*(1:n)+d1[3]*(1:n)^2
    d2 <- c(sample(c(0, runif(1, -3,3)), 1), sample(c(0, runif(1, 0,1)),1))
    X <- u1+T1
    Y <- cumsum(d2[1]+d2[2]*(1:n)+u2)
    Y <- ts(Y/sd(Y)*10, start=start, frequency=frequency)
    X <- ts(X/sd(X)*10, start=start, frequency=frequency)
    sel <- sample(c(0,1), 1)
    if (sel==1) dat <- cbind(X=X,Y=Y) else dat <- cbind(X=Y, Y=X)
    attr(dat, "which") <- c(paste("X is a ", ifelse(sel==1, "TS", "Unit Root")),
                            paste("Y is a ", ifelse(sel==0, "TS", "Unit Root")))
    dat
}

## Demand and Supply

genDemSup <- function(n, a=c(20, 1, -2), d=c(100,-1,2), sige=16,
                    sigeta=16, sigX=5, sigZ=5)
{
    e <- rnorm(n, 0, sqrt(sige))
    eta <- rnorm(n, 0, sqrt(sigeta))
    X <- rnorm(n, 10, sqrt(sigX))
    Z <- rnorm(n, 10, sqrt(sigZ)) 
    Q <- (a[2]*(d[1]+eta+d[3]*X)-d[2]*(a[1]+a[3]*Z+e))/(a[2]-d[2])
    P <- ((d[1]+eta+d[3]*X)-(a[1]+a[3]*Z+e))/(a[2]-d[2])
    dat <- data.frame(P,Q,Z,X,e,eta)
    obj <- list(dat=dat, a=a, d=d, sig=c(e=sige, eta=sigeta),
                sigExo=c(X=sigX, Z=sigZ))
    class(obj) <- "DemSup"
    obj
}
plot.DemSup <- function(x, y=NULL, nCurves=NULL, nPoints=NULL, fac=c(0.9,1.1), ...)
{
    dat <- x$dat
    if (!is.null(nPoints))
        dat <- dat[1:nPoints,]
    plot(Q~P, dat, pch=21, col="lightblue", bg="lightblue",
         main="Demands and Supplies with their equilibrium points",
         bty='n')
    if (is.null(nCurves))
    {
        n <- nrow(dat)
    } else {
        n <- nCurves
        if (nCurves>nrow(dat))
            stop("The number of curves cannot exceed the number of points")
    }
    y <- x    
    for (i in 1:n)
    {
        curve(y$d[1]+y$d[2]*x+y$d[3]*dat$X[i]+dat$eta[i],
              dat$P[i]*fac[1], dat$P[i]*fac[2],
              col="green", add=TRUE)
        curve(y$a[1]+y$a[2]*x+y$a[3]*dat$Z[i]+dat$e[i],
              dat$P[i]*fac[1], dat$P[i]*fac[2],
              col="orange", add=TRUE)
    }
    points(dat$P[1:n], dat$Q[1:n], pch=23, col="darkred",bg="darkred")
    grid()
    invisible()
}

print.DemSup <- function(x, ...)
{
    cat("\\begin{eqnarray*}\n")
    cat("Q^d &=& ", x$d[1], ifelse(x$d[2]<0, "-", "+"), abs(x$d[2]), "P",
        ifelse(x$d[3]<0, "-", "+"), abs(x$d[3]), "X + \\eta\\\\\n",
        "Q^s &=& ", x$a[1], ifelse(x$a[2]<0, "-", "+"), abs(x$a[2]), "P",
        ifelse(x$a[3]<0, "-", "+"), abs(x$a[3]), "Z + e\n",
        "\\end{eqnarray*}\n", sep="")
}

