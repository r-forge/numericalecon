## Extract for tsls

setMethod("extract", "gmmfit", 
          function(model, includeJTest=TRUE, includeFTest=TRUE, ...)
              {
                  s <- summary(model, ...)
                  spec <- modelDims(model@model)
                  coefs <- s@coef
                  names <- rownames(coefs)
                  coef <- coefs[, 1]
                  se <- coefs[, 2]
                  pval <- coefs[, 4]
                  n <- model@model@n
                  gof <- numeric()
                  gof.names <- character()
                  gof.decimal <- logical()
                  if (includeJTest) {
                      if (spec$k == spec$q)
                          {
                              obj.fcn <- NA
                              obj.pv <- NA
                          } else {
                              obj.fcn <- s@specTest@test[1]
                              obj.pv <- s@specTest@test[3]
                          }
                      gof <- c(gof, obj.fcn, obj.pv)
                      gof.names <- c(gof.names, "J-test Statistics", "J-test p-value")
                      gof.decimal <- c(gof.decimal, TRUE, TRUE)
                  }
                  if (includeFTest) {
                      str <- s@strength$strength
                      if (is.null(str))
                          {
                              gof <- c(gof, NA)
                              gof.names <- c(gof.names, "First Stage F-stats")
                              gof.decimal <- c(gof.decimal, TRUE)
                          } else {
                              for (i in 1:nrow(str))
                                  {
                                      gof <- c(gof, str[i,1])
                                      gofn <- paste("First Stage F-stats(",
                                                    rownames(str)[i], ")", sep="")
                                      gof.names <- c(gof.names, gofn)
                                      gof.decimal <- c(gof.decimal, TRUE)
                                  }
                          }
                  }
                  tr <- createTexreg(coef.names = names, coef = coef, se = se, 
                                     pvalues = pval, gof.names = gof.names, gof = gof, 
                                     gof.decimal = gof.decimal)
                  return(tr)
              })

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

genDemSup <- function(n, a=c(20, 1, -2), d=c(100,-1,2,0), sige=16,
                    sigeta=16, sigX=c(5,5), sigZ=5)
{
    e <- rnorm(n, 0, sqrt(sige))
    eta <- rnorm(n, 0, sqrt(sigeta))
    Z <- rnorm(n, 10, sqrt(sigZ))     
    if (d[4] == 0)
    {
        X <- rnorm(n, 10, sqrt(sigX[1]))
        Q <- (a[2]*(d[1]+eta+d[3]*X)-d[2]*(a[1]+a[3]*Z+e))/(a[2]-d[2])
        P <- ((d[1]+eta+d[3]*X)-(a[1]+a[3]*Z+e))/(a[2]-d[2])
        dat <- data.frame(P,Q,Z,X,e,eta)
    } else {
        X1 <- rnorm(n, 10, sqrt(sigX[1]))
        X2 <- rnorm(n, 10, sqrt(sigX[2]))
        Q <- (a[2]*(d[1]+eta+d[3]*X1+d[4]*X2)-d[2]*(a[1]+a[3]*Z+e))/(a[2]-d[2])
        P <- ((d[1]+eta+d[3]*X1+d[4]*X2)-(a[1]+a[3]*Z+e))/(a[2]-d[2])
        dat <- data.frame(P,Q,Z,X1,X2,e,eta)
    }
    obj <- list(dat=dat, a=a, d=d, sig=c(e=sige, eta=sigeta),
                sigExo=list(X=sigX, Z=sigZ))
    class(obj) <- "DemSup"
    obj
}
plot.DemSup <- function(x, y=NULL, nCurves=NULL, nPoints=NULL, fac=c(0.9,1.1), ...)
{
    dat <- x$dat
    if (!is.null(nPoints))
    {
        n <- nPoints
        if (n<2)
            stop("The minimum number of points is 2")
    } else {
        n <- nrow(dat)
    }
    plot(Q~P, dat[1:n,], pch=21, col="lightblue", bg="lightblue",
         main="Demands and Supplies with their equilibrium points",
         bty='n')
    if (!is.null(nCurves))
    {
        if (nCurves < 1)
            stop("nCurves must be greater than 0")
        n <- min(n, nCurves)
    }
    y <- x    
    for (i in 1:n)
    {
        if (x$d[4] == 0)
        {
            curve(y$d[1]+y$d[2]*x+y$d[3]*dat$X[i]+dat$eta[i],
                  dat$P[i]*fac[1], dat$P[i]*fac[2],
                  col="green", add=TRUE)
            curve(y$a[1]+y$a[2]*x+y$a[3]*dat$Z[i]+dat$e[i],
                  dat$P[i]*fac[1], dat$P[i]*fac[2],
                  col="orange", add=TRUE)
        } else {
            curve(y$d[1]+y$d[2]*x+y$d[3]*dat$X1[i]+y$d[4]*dat$X2[i]+dat$eta[i],
                  dat$P[i]*fac[1], dat$P[i]*fac[2],
                  col="green", add=TRUE)
            curve(y$a[1]+y$a[2]*x+y$a[3]*dat$Z[i]+dat$e[i],
                  dat$P[i]*fac[1], dat$P[i]*fac[2],
                  col="orange", add=TRUE)
        }
    }
    points(dat$P[1:n], dat$Q[1:n], pch=23, col="darkred",bg="darkred")
    grid()
    invisible()
}

print.DemSup <- function(x, ...)
{
    cat("\\begin{eqnarray*}\n")
    if (x$d[4] == 0)
        cat("Q^d &=& ", x$d[1], ifelse(x$d[2]<0, "-", "+"), abs(x$d[2]), "P",
            ifelse(x$d[3]<0, "-", "+"), abs(x$d[3]), "X + \\eta\\\\\n", sep="")
    else
        cat("Q^d &=& ", x$d[1], ifelse(x$d[2]<0, "-", "+"), abs(x$d[2]), "P",
            ifelse(x$d[3]<0, "-", "+"), abs(x$d[3]), "X_1",
            ifelse(x$d[4]<0, "-", "+"), abs(x$d[4]), "X_2+\\eta\\\\\n", sep="")
    cat("Q^s &=& ", x$a[1], ifelse(x$a[2]<0, "-", "+"), abs(x$a[2]), "P",
        ifelse(x$a[3]<0, "-", "+"), abs(x$a[3]), "Z + e", sep="")
    cat("\n\\end{eqnarray*}\n")
}

## Computing probabilities from printed Normal tables

probFromTable <- function(a, b=NULL, type=c("less","greater"), digits=3,
                          mu=0, sig2=1, varName="X", approx=FALSE, print=TRUE,
                          CDF=FALSE)
{
    type <- match.arg(type)
    approx <- ifelse(approx, "\\approx", "=")
    if (toupper(varName) == "Z")
        varName <- "X"
    if (is.null(b))
    {
        z <- (a-mu)/sqrt(sig2)
        approxz <- ifelse(z != round(z,2), "\\approx", "=")
        z <- round(z,2)
        ineq <- ifelse(type=="less", "\\leq", "\\geq")
        m <- paste(varName, ineq, a, ")&", approx, "&",
                   "\\Pr\\left(Z", ineq, "\\frac{", a, "-", ifelse(mu<0, "(", ""), mu,
                   ifelse(mu<0, ")", ""),
                   "}{\\sqrt{", sig2, "}}\\right)\\\\\n&", approxz, "&",
                   "Pr(Z", ineq, z, ")\\\\\n")
        if (type == "less")
        {
            if (z > 0)
            {
                if (CDF)
                {
                    m1 <- paste("&\\approx&", round(pnorm(z),digits))
                } else {
                    m1 <- paste("&=&+Pr(0\\leq Z \\leq", z, ")\\\\\n&\\approx&",
                                "0.5+",round(pnorm(z)-.5,digits), "\\\\\n&=&", 0.5+round(pnorm(z)-.5,digits))
                }
                p <- 0.5+round(pnorm(z)-.5,digits)
            } else {
                if (CDF)
                {
                    m1 <- paste("&=&Pr(Z\\geq", -z, ")\\\\\n&=&1-Pr(Z\\leq", -z, ")\\\\\n&\\approx&",
                                "1-", round(pnorm(-z), digits), "\\\\\n&=&", 1-round(pnorm(-z),digits))
                } else {
                    m1 <- paste("&=&0.5-Pr(", z, "\\leq Z \\leq 0)\\\\\n&=&0.5-Pr(0\\leq Z\\leq",
                                -z, ")\\\\\n",
                                "&\\approx&0.5-", round(pnorm(-z)-.5,digits), "\\\\\n&=&",
                                1-round(pnorm(-z),digits))
                }
                p <- 1-round(pnorm(-z),digits)
            }
        } else {
            if (z > 0)
            {
                if (CDF)
                {
                    m1 <- paste("&=&1-Pr(Z\\leq", z, ")\\\\\n&\\approx&",
                                "1-", round(pnorm(z), digits), "\\\\\n&=&", 1-round(pnorm(z),digits))
                } else {                
                    m1 <- paste("&=&0.5-Pr(0\\leq Z \\leq", z, ")\\\\\n&\\approx&",
                                "0.5-",round(pnorm(z)-.5,digits), "\\\\\n&=&", 0.5-round(pnorm(z)-.5,digits))
                }
                p <- 0.5-round(pnorm(z)-.5,digits)
            } else {
                if (CDF)
                {
                    m1 <- paste("&=&1-Pr(Z\\leq", -z, ")\\\\\n&\\approx&",
                                "1-", round(pnorm(-z), digits), "\\\\\n&=&", 1-round(pnorm(-z),digits))
                } else {                                
                    m1 <- paste("&=&0.5+Pr(", z, "\\leq Z \\leq 0)\\\\\n&=&0.5+Pr(0\\leq Z\\leq",
                                -z, ")\\\\\n",
                                "&\\approx&0.5+", round(pnorm(-z)-.5,digits), "\\\\\n&=&",
                                round(pnorm(-z),digits))
                }
                p <- round(pnorm(-z),digits)
            }            
        }
    } else {
        a <- sort(c(a,b))
        z1 <- (a[1]-mu)/sqrt(sig2)
        z2 <- (a[2]-mu)/sqrt(sig2)
        approxz <- ifelse(z1!=round(z1,2) | z2!=round(z2,2), "\\approx", "=")
        z1 <- round(z1,2)
        z2 <- round(z2,2)
        m <- paste(a[1], "\\leq", varName, "\\leq", a[2], ")&", approx, "&",
                   "Pr\\left(\\frac{", a[1], "-", ifelse(mu<0, "(", ""), mu,
                   ifelse(mu<0, ")", ""), "}{\\sqrt{", sig2, "}}",
                   "\\leq Z \\leq",
                   "\\frac{", a[2], "-", ifelse(mu<0, "(", ""), mu, ifelse(mu<0, ")", ""),
                   "}{\\sqrt{", sig2, "}}\\right)",
                   "\\\\\n&", approxz, "&Pr(", z1, "\\leq Z \\leq", z2, ")\\\\\n")
        if (z1==z2)
        {
            m1 <- "&=& 0"
            p <- 0
        } else  if (z1==0) {
            if (CDF)
            {
                m1 <- paste("&=&Pr(Z\\leq", z2,")-0.5\\\\\n&\\approx&",
                            round(pnorm(z2), digits), "-0.5\\\\\n&=&",
                            round(pnorm(z2),digits)-0.5)
            } else {
                m1 <- paste("&\\approx&",round(pnorm(z2)-0.5,digits))
            }
            p <- round(pnorm(z2)-0.5,digits)
        } else if (z2==0) {
            if (CDF)
            {
                m1 <- paste("&=&Pr(0\\leq Z\\leq", -z1, ")\\\\\n",
                            "&=&Pr(Z\\leq", -z1,")-0.5\\\\\n&\\approx&",
                            round(pnorm(-z1), digits), "-0.5\\\\\n&=&",
                            round(pnorm(-z1),digits)-0.5)
            } else {            
                m1 <- paste("&=&Pr(0\\leq Z\\leq", -z1, ")\\\\\n&\\approx&",round(pnorm(z1)-0.5,digits))
            }
            p <- round(pnorm(z1)-0.5, digits)
        } else if (z1<0 & z2<0){
            if (CDF)
            {
                m1 <- paste("&=&Pr(", -z2, "\\leq Z \\leq", -z1, ")\\\\\n",
                            "&=&Pr(Z\\leq", -z1,")-","Pr(Z\\leq", -z2,")\\\\\n&\\approx&",
                            round(pnorm(-z1), digits), "-", round(pnorm(-z2), digits),"\\\\\n&=&",
                            round(pnorm(-z1),digits)-round(pnorm(-z2),digits))
            } else {            
                m1 <- paste("&=&Pr(", z1, "\\leq Z \\leq 0)-Pr(", z2, "\\leq Z \\leq 0)",
                            "\\\\\n&=&",
                            "Pr(0\\leq Z\\leq ", -z1, ")-Pr(0\\leq Z\\leq", -z2, ")",
                            "\\\\\n&\\approx&",round(pnorm(-z1)-0.5,digits), "-", round(pnorm(-z2)-.5,digits),
                            "\\\\\n&=&",round(pnorm(-z1)-0.5,digits)-round(pnorm(-z2)-.5,digits))
            }
            p <- round(pnorm(-z1)-0.5,digits)-round(pnorm(-z2)-.5,digits)
        } else if (z1<0 & z2>0) {
            if (CDF)
            {
                m1 <- paste("&=&Pr(Z\\leq", z2, ")-Pr(Z\\leq", z1, ")\\\\\n",
                            "&=&Pr(Z\\leq", z2, ")-Pr(Z\\geq", -z1, ")\\\\\n",
                            "&=&Pr(Z\\leq", z2, ")-[1-Pr(Z\\leq", -z1, ")]\\\\\n",
                            "\\\\\n&\\approx&",
                            round(pnorm(z2), digits), "-[1-", round(pnorm(-z1), digits),"]\\\\\n&=&",
                            round(pnorm(z2),digits)-(1-round(pnorm(-z1),digits)))
            } else {                        
                m1 <- paste("&=&Pr(", z1, "\\leq Z \\leq 0)+Pr(0\\leq Z\\leq", z2, ")",
                            "\\\\\n&=&",
                            "Pr(0\\leq Z\\leq", -z1,")+Pr(0\\leq Z\\leq", z2, ")",
                            "\\\\\n&\\approx&",round(pnorm(-z1)-0.5,digits), "+", round(pnorm(z2)-.5,digits),
                            "\\\\\n&=&",round(pnorm(-z1)-0.5,digits)+round(pnorm(z2)-.5,digits))
            }
            p <- round(pnorm(-z1)-0.5,digits)+round(pnorm(z2)-.5,digits)
        } else if (z1>0 & z2>0) {
            if (CDF)
            {
                m1 <- paste("&=&Pr(Z\\leq", z2, ")-Pr(Z\\leq", z1, ")\\\\\n",
                            "&\\approx&",
                            round(pnorm(z2), digits), "-", round(pnorm(z1), digits),"\\\\\n&=&",
                            round(pnorm(z2),digits)-round(pnorm(z1),digits))
            } else {                                    
                m1 <- paste("&=&Pr(0\\leq Z\\leq ", z2, ")-Pr(0\\leq Z\\leq", z1, ")",
                            "\\\\\n&\\approx&",
                            round(pnorm(z2)-0.5,digits), "-", round(pnorm(z1)-.5,digits),
                            "\\\\\n&=&",round(pnorm(z2)-0.5,digits)-round(pnorm(z1)-.5,digits))
            }
            p <- round(pnorm(z2)-0.5,digits)-round(pnorm(z1)-.5,digits)
        }        
        
    }
    m <- paste("\\begin{eqnarray*}\nPr(", m, m1, "\n\\end{eqnarray*}\n")
    if (print)
    {
        cat(m)
        return(invisible())
    } else {p}
}


## Generating statistics tables

.chiTab <- function(df=NULL, digits=4, label=NULL)
{
    if(is.null(df))
        df <- seq(1,30)
    if (length(df)<5)
        stop("the number of degrees of freedom must be greater than 4")
    if (is.null(label))
        label <- "chisqTable"
    pr <- c(.01,.025, .05,.1, .9, .95, .975, .99)
    tab <- sapply(pr, function(p) qchisq(p, df))
    tab <- formatC(tab, digits=digits, format="f")
    tab <- cbind(as.character(df),tab)
    colnames(tab) = c("df",formatC(pr,format="f", digits=digits))
    nc <- ncol(tab)
    nr <- nrow(tab)
    com1 <- paste("Example: $x$ that solves $P(Q\\leq x)=0.95$, where $Q$ is ",
                  "$\\chi^2_", df[2], "$, is $x=",
                  round(qchisq(0.95, df[2]),2), "$", sep="")
    com2 <- paste("Example: $x$ that solves $P(Q\\leq x)=0.10$, where $Q$ is ",
                  "$\\chi^2_", df[4], "$, is $x=",
                  round(qchisq(0.10, df[4]),2), "$", sep="")                                 
    com1 <- paste0("\\hline \n \\multicolumn{",nc,"}{l}",
                   "{",com1,"} \n")
    com2 <- paste0("\\\\ \n \\multicolumn{",nc,"}{l}",
                   "{",com2,"} \n")
    addr1 <- list(pos=list(nr), command=paste0(com1,com2))
    t1 <- xtable(tab,caption="Quantiles of the Chi-Square Distribution",
                 label=label)
    list(table=t1, addToRow=addr1)
}

.fTab <- function(df1=NULL, df2=NULL, size=0.05, digits=4, label=NULL)
{
    if (is.null(df1))
        df1 <- 1:10
    if (is.null(df2))
        df2 <- c(seq(10,30),40,60,90,120,Inf)
    if (length(df1)<5 | length(df2)<5)
        stop("The number of df1 and df2 must exceed 4")
    if (is.null(label))
        label <- "fTable"    
    q <- 1-size
    tab <- sapply(df1, function(df) qf(q, df, df2))
    tab <- formatC(tab, format="f", digits=digits)
    tab <- cbind(as.character(df2),tab)
    colnames(tab) = c("df2\\df1",as.character(df1))
    nc <- ncol(tab)
    nr <- nrow(tab)
    com1 <- paste("Example: The ", q*100, "\\% quantile of $F$",
                  "where $F$ is $F(", df1[1], ",", df2[3], ")$ ",
                  "is ", round(qf(q, df1[1], df2[3]),2), sep="")
    com2 <- paste("Example: The ", q*100, "\\% quantile of $F$",
                  "where $F$ is $F(", df1[3], ",", df2[4], ")$ ",
                  "is ", round(qf(q, df1[3], df2[4]),2), sep="")    
    com1 <- paste0("\\hline \n \\multicolumn{",nc,"}{l}",
                   "{",com1,"} \n")
    com2 <- paste0("\\\\ \n \\multicolumn{",nc,"}{l}",
                   "{",com2,"} \n")
    addr1 <- list(pos=list(nr), command=paste0(com1,com2))
    t1 <- xtable(tab,caption="95\\% Quantiles of the F-Distribution",
                 label=label)
    list(table=t1, addToRow=addr1)
}

.tTab <- function(df=NULL, digits=4, label=NULL)
{
    if (is.null(df))
        df <- c(seq(1,30),40,60,90,120, Inf)
    if (length(df)<5)
        stop("the number of degrees of freedom must be greater than 4")
    if (is.null(label))
        label <- "tTable"    
    pr <- c(.90,.95, .975, 0.99, .995)
    tab <- sapply(pr, function(p) qt(p, df))
    tab <- formatC(tab, digits=digits, format="f")
    tab <- cbind(as.character(df),tab)
    colnames(tab) = c("df",formatC(pr,format="f", digits=3))
    nc <- ncol(tab)
    nr <- nrow(tab)
    com1 <- paste("Example: $P(t\\leq x)=0.975$, where $t$ is $t_{", df[2],"}$",
                  " implies $x=$", round(qt(.975, df[2]), digits), sep="")
    com2 <- paste("Example: $P(t\\leq x)=0.99$, where $t$ is $t_{", df[4],"}$",
                  " implies $x=$", round(qt(.99, df[4]), digits), sep="")
    com3 <- paste("Example: $P(t\\leq x)=0.01$, where $t$ is $t_{", df[4],"}$",
                  " implies $x=$", round(qt(.01, df[4]), digits), sep="")        
    
    com1 <- paste0("\\hline \n \\multicolumn{",nc,"}{l}",
                   "{",com1,"} \n")
    com2 <- paste0("\\\\ \n \\multicolumn{",nc,"}{l}",
                   "{",com2,"} \n")
    com3 <- paste0("\\\\ \n \\multicolumn{",nc,"}{l}",
                   "{",com3,"} \n")
    addr1 <- list(pos=list(nr), command=paste0(com1,com2,com3))
    t1 <- xtable(tab,caption="Quantiles of the t-Distribution",
                 label=label)
    list(table=t1, addToRow=addr1)
}

.normTab <- function(digits=4, CDF=FALSE, label=NULL)
{                    
    z <- seq(0,3,by=0.1)
    add <- ifelse(CDF, 0, -0.5)
    tab <- sapply(0:9, function(d) pnorm(z+d/100)+add)
    tab <- formatC(tab, digits=digits, format="f")
    tab <- cbind(formatC(z,digits=1, format="f"),tab)
    colnames(tab) = c("z",seq(0,9))
    nc <- ncol(tab)
    nr <- nrow(tab)
    if (is.null(label))
        label <- "normalTable"    
    com1 <- "Example 1: If $Z\\sim N(0,1)$, then $P(0\\leq Z\\leq 1.84)=0.467$"
    com2 <- "Example 2:  $P(-1.54 \\leq Z \\leq 0)=P(0\\leq Z\\leq 1.54) = 0.438$"
    com1 <- paste0("\\hline \n \\multicolumn{",nc,"}{l}",
                   "{",com1,"} \n")
    com2 <- paste0("\\\\ \n \\multicolumn{",nc,"}{l}",
                   "{",com2,"} \n")
    addr1 <- list(pos=list(nr), command=paste0(com1,com2))
    t1 <- xtable(tab,caption="$Pr(0\\leq Z\\leq z)$, where $Z\\sim N(0,1)$",
                 label=label)
    list(table=t1, addToRow=addr1)
}


distTable <- function(df=NULL, df2=NULL, digits=4, dist=c("t", "chisq", "f", "normal"),
                      size=0.05, CDF=FALSE, float=TRUE, label=NULL)
{
    dist <- match.arg(dist)
    tab <- switch(dist,
                  t = get(".tTab")(df, digits, label),
                  chisq = get(".chiTab")(df, digits, label),
                  f = get(".fTab")(df, df2, size, digits, label),
                  normal = get(".normTab")(digits, CDF, label))
    print(tab$table, include.rownames=FALSE, caption.placement="top",
          hline.after=c(-1, 0), add.to.row = tab$addToRow, comment=FALSE,
          floating=float)
}
