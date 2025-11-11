## Print Tools

printReg <- function(obj,digits=4, format='f', form=NULL, maxpl=5, adjrsq=FALSE, 
                     se=NULL, stars=FALSE, dist=c("t","n"), label=NULL,
                     omit=NULL, ...)
    {
        dist <- match.arg(dist)
        if (is.null(label))
            {
                typeEQ <- "equation*"
                label=""
            } else {
                typeEQ <- "equation"
                label <- paste("\\label{",label,"}", sep="")
            }
        if (!is.null(omit))
            {
                omit <- omit[omit != "(Intercept)"]
                chk <- lapply(omit, function(o) grep(o, names(obj$coef)))
                chk <- sapply(1:length(chk), function(i) length(chk[[i]])==0)
                omit <- omit[!chk]
                if (length(omit)==0)
                    omit <- NULL
            }
        if (is.null(form))
            {                
                cat("\\begin{",typeEQ, "}", label,"\n", sep="")
                cat("\\begin{split}\n")
                ncoef <- names(coef(obj))
                if (!is.null(omit))
                    {
                        omit <- do.call("c", lapply(omit, function(o) grep(o, ncoef)))
                        omit <- -unique(omit)
                        nomit <- length(omit)
                        add <- paste("\\mbox{ (+ ", nomit, " omitted terms)}", sep="")
                    } else {
                        add <- ""
                        omit <- 1:length(ncoef)
                    }                
                Intercept <- attr(obj$terms, "intercept")
                b <- formatC(abs(obj$coef), digits=digits, format=format, ...)[omit]
                if (!is.null(se))
                    snum <- se[omit]
                else
                    snum <- summary(obj)$coef[omit,2]
                s <- formatC(snum, digits=digits, format=format, ...)
                ncoef <- ncoef[omit]
                if (stars)
                    {
                        ttest <- coef(obj)[omit]/snum
                        if (dist=="t")
                            pv <- 2*pt(-abs(ttest), obj$df.residual)
                        else
                            pv <- 2*pnorm(-abs(ttest))
                        sym <- symnum(pv, cutpoints=c(0,.01,.05,.1,1),
                                      symbols=c("^{***}","^{**}","^*"," "))
                        symmess <- "\\\\& ^*\\text{pv}<0.1\\mbox{; }^{**}\\text{pv}<0.05\\mbox{; }^{***}\\text{pv}<0.01"
                    } else {
                        sym <- rep("", length(coef(obj)[omit]))
                        symmess <- ""
                    }
                if (length(attr(obj$terms, "factors")))
                {
                    ny <- rownames(attr(obj$terms, "factors"))[1]                     
                } else {
                    ny <- as.character(formula(obj))[2]
                }
                ny <- paste("\\widehat{",ny,"}",sep="")
                if (inherits(obj, "glm"))
                    ny <- paste("\\text{link}\\left[", ny, "\\right]", sep="")
                cat(ny,"&=")
                if (obj$coef[1] < 0)
                    cat("\\underset{(",s[1],")",sym[[1]],"}{-",b[1],"}", sep="")
                else
                    cat("\\underset{(",s[1],")",sym[[1]],"}{",b[1],"}")                    
                if (Intercept==0)
                        cat("~",ncoef[1], sep="")
                j <- 1
                if (length(b)>=2)
                {
                    for (i in 2:length(b))
                    {
                        if (j>maxpl)
                        {
                            j <- 1
                            cat("\\\\&\\quad\n")
                        }
                        if ((obj$coef[omit])[i] < 0)
                            cat("~-~")
                        else
                            cat("~+~")
                        cat("\\underset{(",s[i],")",sym[[i]],"}{",b[i],"}~")
                        cat(ncoef[i])
                        j <- j+1
                    }
                }
                cat(add)
                n <- length(obj$residuals)
                if (inherits(obj, "glm"))
                {
                    cat("\\\\ &\\quad n=", n, ",~~AIC=", round(obj$aic,digits))
                    cat(",~~\\mbox{Residual Deviance}=", round(obj$deviance,digits))
                    cat(",\\\\ &\\quad \\mbox{Null Deviance}=",
                        round(obj$null.deviance,digits))
                    cat(",~~\\mbox{Family}=\\mbox{", obj$family$family,
                        "},~~\\mbox{Link}=\\mbox{", obj$family$link, "}", sep="")
                } else {
                    cat("\\\\ &\\quad n=", n, ",~~R^2=", round(summary(obj)$r.squared,digits))
                    cat(", SSR=", round(sum(obj$resid^2),digits))
                    if (adjrsq)
                        cat(", \\bar{R}^2=", round(summary(obj)$adj,digits))
                }
                if (!is.null(se))
                    cat("\\mbox{ (Robust S-E)}")
                cat(symmess)
                cat("\n\\end{split}\n")
                cat("\\end{", typeEQ, "}\n", sep="")
            } else {
                t <- terms(form)
                y <- rownames(attr(t, "factors"))[1]                     
                x <- colnames(attr(t, "factors"))
                cat("\\begin{",typeEQ, "}", label,"\n\\begin{split}\n",
                    y,"&=\\beta_0", sep="")
                j  <-  1
                for (i in 1:length(x))
                    {
                        if (j>maxpl)
                            {
                                j <- 1
                                cat("\\\\&\n")
                            }
                        cat("+\\beta_",i,x[i],sep="")
                        j <- j+1
                    }
                cat("+u\n\\end{split}\n\\end{",typeEQ,"}", sep="")
            }
    }

## print method for solution generators

print.metricsSol <- function(x, addMess=NULL, ...)
{
    cat("\n",x$test,"\n\n")
    cat("\n",x$ftest,"\n")
    if (!is.null(addMess))
        x$mes <- paste(x$mes, " \\textbf{(", addMess, ")}", sep="")
    cat(x$mes, "\n")
}


## Solution generators: Inference on the means

testMean <- function(x, h0, size=.05, alter=c("diff","greater","less"),
                     assume=c("Normal","nonNormal"), digits=4, dround=NULL,
                     xbar=NULL, se=NULL, n=NULL)
{    
    alter <- match.arg(alter)
    assume <- match.arg(assume)
    if (!is.null(xbar))
    {
        if (is.null(se) | is.null(n))
            stop("the standard error and/or sample size is missing")
        m <- xbar
        s = se
    } else {
        m <- mean(x)
        s <- sd(x)
        if (!is.null(dround))
        {
            m <- round(m, dround)
            s <- round(s, dround)
        }
        n <- length(x)
    }
    what <- paste("Testing $H_0:~\\mu=",h0,"$ against $H_1:~\\mu",
                  switch(alter, diff="\\neq", greater=">", less="<"),h0, "$",
                  " at ", size*100, "\\%", sep="")    
    t <- (m-h0)/(s/sqrt(n))  
    if (assume=="Normal")
    {
        dist <- paste("t_{",n-1,"}",sep="")
        mes2 <- ""
        distAs <- "\\sim"
    } else {
        dist <- "N(0,1)"
        mes2 <- " (The N(0,1) is an approximation based on the C.L.T. because the distribution of the data is unknown)"
        distAs <- "\\approx"
    }
    ftest <- paste("\\[\ntest = \\frac{(\\bar{x}-",h0,")}{s/\\sqrt{n}} = \\frac{(",
                   round(m,digits),"-",h0,")}{(",round(s,digits),")/\\sqrt{",
                   n,"}}=",round(t,digits),distAs, dist,"\n\\]")               
    if (alter == "diff")
    {
        crit <- ifelse(assume=="nonNormal", qnorm(1-size/2), qt(1-size/2,n-1))
        if (abs(t)>crit)
            mes <- paste("Since $|", round(t,digits), "|>", round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $|", round(t,digits), "|<", round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else if (alter == "greater") {
        crit <- ifelse(assume=="nonNormal", qnorm(1-size), qt(1-size,n-1))
        if (t>crit)
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else {
        crit <- ifelse(assume=="nonNormal", qnorm(size), qt(size,n-1))
        if (t<crit)
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    }
    mes <- paste(mes, mes2, sep="")
    obj <- list(ftest=ftest, mes=mes, type="Test on the mean", test=what)
    class(obj) <- "metricsSol"
    obj
}        

ciMean <- function(x, size=0.05, assume=c("Normal","nonNormal"), digits=4, dround=NULL,
                   xbar=NULL, se=NULL, n=NULL)
{
    assume <- match.arg(assume)
    if (!is.null(xbar))
    {
        if (is.null(se) & is.null(n))
            stop("the standard error and/or sample size is missing")
        m <- xbar
        s = se
    } else {
        m <- mean(x)
        s <- sd(x)
        if (!is.null(dround))
        {
            m <- round(m, dround)
            s <- round(s, dround)
        }
        n <- length(x)
    }
    what <- paste((1-size)*100, "\\% confidence interval for the mean", sep="")        
    crit <- ifelse(assume=="nonNormal", qnorm(1-size/2), qt(1-size/2,n-1))
    up <- m+crit*s/sqrt(n)
    down <- m-crit*s/sqrt(n)
    ftest <- paste("\\[\n \\begin{split}\n",
                   "CI & = \\left[ \\bar{x}-t^*\\frac{s}{\\sqrt{n}}",
                   ",~\\bar{x}+t^*\\frac{s}{\\sqrt{n}}\\right]",
                   "\\\\ \n& = \\left[", round(m,digits),"-",round(crit,digits),
                   "\\frac{",round(s,digits),"}{\\sqrt{",n,"}}",
                   ",~", round(m,digits),"+",round(crit,digits),
                   "\\frac{",round(s,digits),"}{\\sqrt{",n,"}}\\right]",
                   "\\\\ \n& = [", round(down,digits),",~",round(up,digits),"]",
                   "\\end{split}\n\\]",sep="")
    if (assume=="Normal")
        mes <- paste("The $t^*$ is the ", (1-size/2)*100, "\\% quantile of the ",
                     "t-distribution with ", n-1, " degrees of freedom.", sep="")
    else
        mes <- paste("The $t^*$ is the ", (1-size/2)*100, "\\% quantile of the ",
                     "$N(0,1)$ (The $N(0,1)$ is an approximation based on the ",
                     "C.L.T. because the distribution of the data is unknown).", 
                     sep="")
    obj <- list(ftest=ftest, mes=mes, type="Confidence interval for the mean",
                test=what)
    class(obj) <- "metricsSol"
    obj
}


testDiffMeans <- function(x1, x2, h0, size=.05, alter=c("diff","greater","less"),
                          assume=c("Normal","nonNormal"), assumev=c("same", "diff"),
                          digits=4, dround=NULL, xbar=NULL, se=NULL, n=NULL)
{
    alter <- match.arg(alter)
    assume <- match.arg(assume)
    assumev <- match.arg(assumev)
    if (!is.null(xbar))
    {
        if (is.null(se) & is.null(n))
            stop("the standard error and/or sample size is missing")
        if (any(c(length(xbar)!=2, length(se)!=2, length(n)!=2)))
            stop("xbar, se and n must contain two numeric values")
        m1 <- xbar[1]; m2 <- xbar[2]
        s1 <- se[1]; s2 <- se[2]
        n1 <- n[1]; n2 <- n[2]
    } else {        
        m1 <- mean(x1)
        m2 <- mean(x2)
        s1 <- sd(x1)
        s2 <- sd(x2)
        if (!is.null(dround))
        {
            m1 <- round(m1, dround)
            m2 <- round(m2, dround)
            s1 <- round(s1, dround)
            s2 <- round(s1, dround)
        }        
        n1 <- length(x1)
        n2 <- length(x2)
    }
    what <- paste("Testing $H_0:~\\mu_1-\\mu_2=",h0,"$ against $H_1:~\\mu_1-\\mu_2",
                  switch(alter, diff="\\neq", greater=">", less="<"),h0, "$",
                  " at ", size*100, "\\%", sep="")    
    if (assumev == "diff")
    {
        s <- sqrt(s1^2/n1+s2^2/n2)
        fs <- paste("\\[s = \\sqrt{\\frac{\\hat{\\sigma}^2_1}{n_1}",
                    "+\\frac{\\hat{\\sigma}^2_2}{n_2}} = ",
                    "\\sqrt{\\frac{",round(s1^2,digits),"}{",n1,"}+",
                    "\\frac{",round(s2^2,digits),"}{",n2,"}}=",round(s,digits),
                    "\n\\]\n",sep="")
    } else {
        s <- sqrt((s1^2*(n1-1)+s2^2*(n2-1))/(n1+n2-2)*(1/n1+1/n2))
        fs <- paste("\\[\n\\begin{split}\n",
                    "s &= \\sqrt{\\frac{1}{n_1+n_2-2}\\left[",
                    "\\hat\\sigma^2_1(n_1-1)",
                    "+ \\hat\\sigma^2_2(n_2-1)\\right]",
                    "\\left(\\frac{1}{n_1}+\\frac{1}{n_2}\\right)}",
                    " \\\\ \n &= \\sqrt{",
                    "\\frac{1}{",n1+n2-2,"}\\left[",round(s1^2,digits), "\\times",
                    "(", n1, "-1) + ", round(s2^2, digits),
                    "\\times (", n2, "-1)\\right]\\left(\\frac{1}{",
                    n1,"}+\\frac{1}{",n2,"}\\right)}",
                    "\\\\ \n &=", round(s,digits=4),"\n\\end{split}",
                    "\n \\]\n",sep="")                        
    }
    t <- (m1-m2-h0)/s
    useT <- assume=="Normal" & assumev=="same"
    if (useT)
    {
        dist <- paste("t_{",n1+n2-2,"}",sep="")
        mes2 <- ""
        distAs <- "\\sim"
    } else {
        dist <- "N(0,1)"
        mes2 <- " (The N(0,1) is an approximation based on the C.L.T. because the distribution of the data is unknown)"
        distAs <- "\\approx"
    }
    rhs <- if(h0==0) {""} else if(h0<0) {paste("+",-h0,sep="")} else {paste("-",h0,sep="")}
    ftest <- paste("\\[\ntest = \\frac{(\\bar{x}_1-\\bar{x}_2)", rhs,
                   "}{s} = \\frac{(",
                   round(m1,digits),"-",round(m2,digits),")", rhs, "}{",round(s,digits),"}=",
                   round(t,digits),distAs, dist,"\n\\]")               
    if (alter == "diff")
    {
        crit <- ifelse(!useT, qnorm(1-size/2), qt(1-size/2,n1+n2-2))
        if (abs(t)>crit)
            mes <- paste("Since $|", round(t,digits), "|>", round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $|", round(t,digits), "|<", round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else if (alter == "greater") {                
        crit <- ifelse(!useT, qnorm(1-size), qt(1-size,n1+n2-2))
        if (t>crit)
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else {
        crit <- ifelse(!useT, qnorm(size), qt(size,n1+n2-2))
        if (t<crit)
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    }
    mes <- paste(mes, mes2, sep="")
    obj <- list(ftest=paste(fs,ftest), mes=mes, type="Test on the difference between means",
                test=what)
    class(obj) <- "metricsSol"
    obj
}        


ciDiffMeans <- function(x1, x2, size=.05, assume=c("Normal","nonNormal"),
                        assumev=c("same", "diff"),
                        digits=4, dround=NULL, xbar=NULL, se=NULL, n=NULL)
{
    assume <- match.arg(assume)
    assumev <- match.arg(assumev)
    if (!is.null(xbar))
    {
        if (is.null(se) & is.null(n))
            stop("the standard error and/or sample size is missing")
        if (any(c(length(xbar)!=2, length(se)!=2, length(n)!=2)))
            stop("xbar, se and n must contain two numeric values")
        m1 <- xbar[1]; m2 <- xbar[2]
        s1 <- se[1]; s2 <- se[2]
        n1 <- n[1]; n2 <- n[2]
    } else {        
        m1 <- mean(x1)
        m2 <- mean(x2)
        s1 <- sd(x1)
        s2 <- sd(x2)
        if (!is.null(dround))
        {
            m1 <- round(m1, dround)
            m2 <- round(m2, dround)
            s1 <- round(s1, dround)
            s2 <- round(s2, dround)
        }        
        n1 <- length(x1)
        n2 <- length(x2)
    }
    what <- paste((1-size)*100, "\\% confidence interval for $(\\mu_1-\\mu_2)$", sep="")        
    if (assumev == "diff")
    {
        s <- sqrt(s1^2/n1+s2^2/n2)
        fs <- paste("\\[s = \\sqrt{\\frac{\\hat{\\sigma}^2_1}{n_1}",
                    "+\\frac{\\hat{\\sigma}^2_2}{n_2}} = ",
                    "\\sqrt{\\frac{",round(s1^2,digits),"}{",n1,"}+",
                    "\\frac{",round(s2^2,digits),"}{",n2,"}}=",round(s,digits),
                    "\n\\]\n",sep="")
    } else {
        s <- sqrt((s1^2*(n1-1)+s2^2*(n2-1))/(n1+n2-2)*(1/n1+1/n2))
        fs <- paste("\\[\n\\begin{split}\n",
                    "s &= \\sqrt{\\frac{1}{n_1+n_2-2}\\left[",
                    "\\hat\\sigma^2_1(n_1-1)",
                    "+ \\hat\\sigma^2_2(n_2-1)\\right]",
                    "\\left(\\frac{1}{n_1}+\\frac{1}{n_2}\\right)}",
                    " \\\\ \n &= \\sqrt{",
                    "\\frac{1}{",n1+n2-2,"}\\left[",round(s1^2,digits), "\\times",
                    "(", n1, "-1) + ", round(s2^2, digits),
                    "\\times (", n2, "-1)\\right]\\left(\\frac{1}{",
                    n1,"}+\\frac{1}{",n2,"}\\right)}",
                    "\\\\ \n &=", round(s,digits=4),"\n\\end{split}",
                    "\n \\]\n",sep="")                        
    }
    useT <- assume=="Normal" & assumev=="same"
    if (useT)
    {
        crit <- qt(1-size/2,n1+n2-2)
        mes <- paste("The $t^*$ is the ", (1-size/2)*100, "\\% quantile of the ",
                     "t-distribution with ", n1+n2-2, " degrees of freedom.", sep="")
    } else {
        crit <- qnorm(1-size/2)
        mes <- paste("The $t^*$ is the ", (1-size/2)*100, "\\% quantile of the ",
                     "$N(0,1)$ (The $N(0,1)$ is an approximation based on the ",
                     "C.L.T. because the distribution of statistic is unknown).", 
                     sep="")
    }
    up <- m1-m2+crit*s
    down <- m1-m2-crit*s
    ftest <- paste("\\[\n \\begin{split}\n",
                   "CI & = \\left[(\\bar{X}_1-\\bar{X}_2)-t^*s",
                   ",~(\\bar{X}_1-\\bar{X}_2)+t^*s\\right]",
                   "\\\\ \n& = \\left[(", round(m1,digits),"-", round(m2,digits), ")-",
                   round(crit,digits), "\\times",round(s,digits),
                   ",~(", round(m1,digits), "-", round(m2,digits), ")+",round(crit,digits),
                   "\\times", round(s,digits),"\\right]",
                   "\\\\ \n& = [", round(down,digits),",~",round(up,digits),"]",
                   "\\end{split}\n\\]",sep="")
    obj <- list(ftest=paste(fs,ftest), mes=mes,
                type="Confidence interval for the difference between means",
                test=what)
    class(obj) <- "metricsSol"
    obj    
}        

## Solution generators: Inference on the variance

testVar <- function(x, h0, size=.05, alter=c("diff","greater","less"),
                    assume=c("Normal","nonNormal"), digits=4, dround=NULL,
                    se=NULL, n=NULL)
{
    alter <- match.arg(alter)
    assume <- match.arg(assume)
    if (!is.null(se))
    {
        if (is.null(n))
            stop("the sample size is missing")
        v <- se^2
    } else {
        v <- var(x)
        if (!is.null(dround))
        {
            v <- round(v, dround)
        }
        n <- length(x)
    }
    t <- (n-1)*v/h0
    what <- paste("Testing $H_0:~\\sigma^2=",h0,"$ against $H_1:~\\sigma^2",
                  switch(alter, diff="\\neq", greater=">", less="<"),h0, "$",
                  " at ", size*100, "\\%", sep="")            
    dist <- paste("\\chi^2_{",n-1,"}",sep="")   
    if (assume=="Normal")
        mes2 <- ""
    else 
        mes2 <- paste(" (The $", dist, 
                      "$ is not valid in this case ",
                      "because the distribution of the data is unknown)", sep="")
    ftest <- paste("\\[\ntest = \\frac{(n-1)\\hat{\\sigma}^2}{",h0,"} = \\frac{(",
                   n-1,")",round(v,digits),"}{",h0,"} = ",
                   round(t,digits),"\\sim", dist,"\n\\]")               
    if (alter == "diff")
    {
        crit1 <- qchisq(size/2, n-1)
        crit2 <- qchisq(1-size/2, n-1)
        if (t>crit2)
            mes <- paste("Since $", round(t,digits), ">", round(crit2,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else if (t<crit1)
            mes <- paste("Since $", round(t,digits), "<", round(crit1,digits),
                         "$ (the ", (size/2)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), ">", round(crit1,digits),
                         "$ (the ", (size/2)*100, "\\% quantile of the $",
                         dist, "$) and $", round(t,digits), "<",
                         round(crit2,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")                
    } else if (alter == "greater") {
        crit <- qchisq(1-size, n-1)
        if (t>crit)
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else {
        crit <- qchisq(size,n-1)
        if (t<crit)
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    }
    mes <- paste(mes, mes2, sep="")
    obj <- list(ftest=ftest, mes=mes,
                type="Test on the variance",
                test=what)
    class(obj) <- "metricsSol"
    obj            
}

ciVar <- function(x, size=0.05, assume=c("Normal","nonNormal"), digits=4, dround=NULL,
                  se=NULL, n=NULL)
{
    assume <- match.arg(assume)    
    if (!is.null(se))
    {
        if (is.null(n))
            stop("the sample size is missing")
        v <- se^2
    } else {
        v <- var(x)
        if (!is.null(dround))
        {
            v <- round(v, dround)
        }
        n <- length(x)
    }
    qu <- qchisq(1-size/2, n-1)
    ql <- qchisq(size/2, n-1)  
    up <- (n-1)*v/ql
    down <- (n-1)*v/qu
    what <- paste((1-size)*100, "\\% confidence interval for the variance", sep="")
    ftest <- paste("\\[\n \\begin{split}\n",
                   "CI & = \\left[\\frac{(n-1)S^2}{Q^*_u}, ",
                   "\\frac{(n-1)S^2}{Q^*_l}\\right]",
                   "\\\\ \n& = \\left[", (n-1),"\\times",
                   "\\frac{",round(v,digits),"}{",round(qu, digits),"}",
                   ",~",(n-1),"\\times",
                   "\\frac{",round(v,digits),"}{",round(ql, digits),"}\\right]",
                   "\\\\ \n& = [", round(down,digits),",~",round(up,digits),"]",
                   "\\end{split}\n\\]",sep="")
    mes <- paste("The $Q^*_u$ is the ", (1-size/2)*100, "\\% quantile of the ",
                 "$\\chi^2_{", n-1, "}$ and $Q^*_l$ is its ", size/2*100,
                 "\\% quantile.", 
                 sep="")
    if (assume=="nonNormal")
        mes <- paste(mes, " (The $\\chi^2_{", n-1, "}$", 
                     "$ is not valid in this case ",
                     "because the distribution of the data is unknown)", sep="")
    obj <- list(ftest=ftest, mes=mes,
                type="Confidence interval for the variance",
                test=what)
    class(obj) <- "metricsSol"
    obj                
}        

testDiffVar <- function(x1, x2, size=.05, 
                        assume=c("Normal","nonNormal"), digits=4, dround=NULL,
                        se=NULL, n=NULL)
{
    assume <- match.arg(assume)
    if (!is.null(se))
    {
        if (is.null(n))
            stop("the sample size is missing")
        if (any(c(length(se)!=2, length(n)!=2)))
            stop("se and n must contain two numeric values")
        v1 <- se[1]^2; v2 <- se[2]^2
        n1 <- n[1]; n2 <- n[2]
    } else {        
        s1 <- sd(x1)
        s2 <- sd(x2)
        if (!is.null(dround))
        {
            s1 <- round(s1, dround)
            s2 <- round(s2, dround)
        }        
        n1 <- length(x1)
        n2 <- length(x2)
        v1 <- s1^2; v2 <- s2^2
    }
    num <- max(v1,v2)
    den <- min(v1,v2)
    t <- num/den
    nnum <- ifelse(v1>v2, n1, n2)
    nden <- ifelse(v1>v2, n2, n1)
    dist <- paste("F(",nnum-1,",",nden-1,")",sep="")
    what <- paste("Testing $H_0:~\\sigma^2_1=\\sigma^2_2$, against $H_1:~\\sigma^2_1",
                  "\\neq \\sigma^2_2$", 
                  " at ", size*100, "\\%", sep="")                
    if (assume=="Normal")
        mes2 <- ""
    else 
        mes2 <- paste(" (The ", dist, 
                      " is not a valid distribution ",
                      "because the distribution of the data is unknown)", sep="")
    ftest <- paste("\\[\ntest = \\frac{\\hat{\\sigma}^2_1}{\\hat{\\sigma}^2_2}",
                   " = \\frac{",round(num,digits),"}{",round(den,digits),"} = ",
                   round(t,digits),"\\sim", dist,"\n\\]")               
    crit <- qf(1-size, nnum-1,nden-1)
    if (t>crit)
        mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                     "$ (the ", (1-size)*100, "\\% quantile of the $",
                     dist, "$), we reject $H_0$.", sep="")
    else
        mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                     "$ (the ", (1-size)*100, "\\% quantile of the $",
                     dist, "$), we do not reject $H_0$.", sep="")
    mes <- paste(mes, mes2, sep="")
    obj <- list(ftest=ftest, mes=mes,
                type="Test on the equality of two variances",
                test=what)
    class(obj) <- "metricsSol"
    obj                
}        

## Solution generators: Inference on regression coefficients

testRegCoef <- function(object, h0, which, size=.05, alter=c("diff","greater","less"),
                        digits=4, assume=c("Normal","nonNormal"),
                        beta=NULL, se = NULL, n=NULL, ncoef=NULL, dround=NULL)
{
    alter <- match.arg(alter)
    assume <-  match.arg(assume)
    signif <- h0==0 & alter=="diff"
    exact <- assume=="Normal" & is.null(beta) & is.null(se)
    if (!is.null(beta))
    {
        if (is.null(se) | is.null(n) | is.null(ncoef))
            stop("the standard error and/or sample size and/or the number of coefficients is missing")
        m <- beta
        s <- se
        df <- n-ncoef
    } else {
        if (missing(which))
            stop("You must specify which coefficient you are performing the test on")
        m <- object$coef[which]
        if (is.null(se))
        {
            s <- sqrt(diag(vcov(object)))[which]
        } else {
            if (length(se)!=length(coef(object)))
                stop("se must be a vector of the same length as the coefficient vector")
            s <- se[which]
        }
        n <- nobs(object)
        df <- object$df
        if (!is.null(dround))
        {
            m <- round(m, dround)
            s <- round(s, dround)
        }
    }
    coefn <- ifelse(!is.null(beta), "\\beta_j", paste("\\beta_", which-1, sep=""))
    what <- paste("Testing $H_0:~", coefn, "=",h0, "$ against $H_1:~", coefn,
                  switch(alter, diff="\\neq", greater=">", less="<"),h0, "$",
                  " at ", size*100, "\\%", sep="")    
    t <- (m-h0)/s
    if (exact)
    {
        dist <- paste("t_{",df,"}",sep="")
        mes2 <- ""
        distAs <- "\\sim"
    } else {
        dist <- "N(0,1)"
        
        mes2 <- " (The N(0,1) is an approximation based on the C.L.T. because either the distribution of the error term is unknown and/or the errors are not homoskedastic.)"
        distAs <- "\\approx"
    }
    ftest <- paste("\\[\ntest = \\frac{(",coefn,"-",h0,")}{se} = \\frac{(",
                   round(m,digits),"-",h0,")}{",round(s,digits),"}=",
                   round(t,digits),distAs, dist,"\n\\]")               
    if (alter == "diff")
    {
        crit <- ifelse(!exact, qnorm(1-size/2), qt(1-size/2,df))
        if (abs(t)>crit)
        {
            mes <- paste("Since $|", round(t,digits), "|>",
                         round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
            if (signif)
                mes <- paste(mes,
                             "The coefficient is therefore significant")
        } else {
            mes <- paste("Since $|", round(t,digits), "|<",
                         round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
            if (signif)
                mes <- paste(mes,
                             "The coefficient is therefore non-significant")
        }
    } else if (alter == "greater") {
        crit <- ifelse(!exact, qnorm(1-size), qt(1-size,df)) 
        if (t>crit)
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else {
        crit <- ifelse(!exact, qnorm(size), qt(size,df))
        if (t<crit)
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    }
    mes <- paste(mes, mes2, sep="")
    obj <- list(ftest=ftest, mes=mes,
                type="Test on one regression coefficient",
                test=what)
    class(obj) <- "metricsSol"
    obj                    
}        


ciRegCoef <- function(object, size=0.05,  which, digits=4,
                      assume=c("Normal", "nonNormal"),
                      beta=NULL, se=NULL, n=NULL, ncoef=NULL, dround=NULL)
{
    assume <- match.arg(assume)
    exact <- assume=="Normal" & is.null(beta) & is.null(se)
    if (!is.null(beta))
    {
        if (is.null(se) | is.null(n) | is.null(ncoef))
            stop("the standard error and/or sample size and/or the number of coefficients is missing")
        m <- beta
        s <- se
        df <- n-ncoef
    } else {
        if (missing(which))
            stop("You must specify for which coefficient you want to compute the confidence interval")
        m <- object$coef[which]
        if (is.null(se))
        {
            s <- sqrt(diag(vcov(object)))[which]
        } else {
            if (length(se)!=length(coef(object)))
                stop("se must be a vector of the same length as the coefficient vector")
            s <- se[which]
        }
        n <- nobs(object)
        df <- object$df
        if (!is.null(dround))
        {
            m <- round(m, dround)
            s <- round(s, dround)
        }
    }
    coefn <- ifelse(!is.null(beta), "\\beta_j", paste("\\beta_", which-1, sep=""))
    what <- paste((1-size)*100, "\\% confidence interval for $", coefn, "$.", sep="")
    crit <- ifelse(exact, qt(1-size/2,df), qnorm(1-size/2))
    dist <- ifelse(exact,
                   paste("t-distribution with ", df, " degrees of freedom.",sep=""),
                   "N(0,1).")
    up <- m+crit*s
    down <- m-crit*s
    ftest <- paste("\\[\n \\begin{split}\n",
                   "CI & = \\left[",coefn,"-t^*se",
                   ",~", coefn, "+t^*se\\right]",
                   "\\\\ \n& = \\left[", round(m,digits),"-",round(crit,digits),
                   "\\times",round(s,digits),
                   ",~", round(m,digits),"+",round(crit,digits),
                   "\\times",round(s,digits),"\\right]",
                   "\\\\ \n& = [", round(down,digits),",~",round(up,digits),"]",
                   "\\end{split}\n\\]",sep="")
    mes <- paste("The $t^*$ is the ", (1-size/2)*100, "\\% quantile of the ",
                 dist, sep="")
    if (!exact)
        mes <- paste(mes, " (The N(0,1) is an approximation based on the C.L.T. ",
                     "because either the distribution of the error term is unknown ",
                     "and/or the errors are not homoskedastic.)", sep="")
    obj <- list(ftest=ftest, mes=mes,
                type="Confidence interval on one regression coefficient",
                test=what)
    class(obj) <- "metricsSol"
    obj                        
}

testLinRegCoefs <- function(object, R, q,  size=.05, alter=c("diff","greater","less"),
                      digits=4, assume=c("Normal","nonNormal"),
                      beta=NULL, vcoef=NULL, n=NULL, dround=NULL)
{        
    alter <- match.arg(alter)
    assume <- match.arg(assume)
    exact <- assume=="Normal" & is.null(beta) & is.null(vcoef)
    R <- c(R)
    q <- c(q)
    if (length(q)!=1)
        stop("The right-hand side must be a scalar")
    if (!is.null(beta))
    {
        if (is.null(vcoef) | is.null(n))
            stop("the covariance matrix and/or sample size is missing")
        b <- beta
        v <- vcoef
        df <- n-length(b)
    } else {
        b <- coef(object)
        if (is.null(vcoef))
        {
            v <- vcov(object)
        } else {
            v <- vcoef
        }
        n <- nobs(object)
        if (!is.null(dround))
        {
            b <- round(b, dround)
            v <- round(v, dround)
        }
        df <- object$df
    }
    if (length(R)!=length(b))
        stop("The number of columns of R does not match the number of coefficients")
    m <- sum(b*R)-q
    s <- sqrt(c(t(R)%*%v%*%R))
    which <- which(R!=0)
    sR <- ifelse(R[which]<0,"-", "+")
    sq <- ifelse(sign(q)==-1, "+", "-")
    if (sR[1] == "+")
        sR[1] <- ""
    mR <- abs(R[which])
    mR[mR==1] <- ""
    coefn <- vector()
    sdevn <- vector()
    for (i in 1:length(mR))
    {
        coefn2 <- paste(sR[i], mR[i], "\\hat{\\beta}_",
                        which[i]-1, sep="")
        coefn <- paste(coefn, coefn2, sep="")
        if (i == 1)
            start <- ""
        else
            start <- "+"
        if (abs(R[which[i]]) == 1)
        {
            sdevn2 <- paste(start, "Var(\\hat{\\beta}_", which[i]-1,")",
                            sep="")
        } else {
            sdevn2 <- paste(start, abs(R[which[i]]), "^2Var(\\hat{\\beta}_",
                            which[i]-1,")",sep="")
        }
        sdevn <- paste(sdevn, sdevn2, sep="")
    }
    coefnW <- gsub("\\\\hat", "", coefn)
    coefn <- paste(coefn, sq, abs(q), sep="") 
    if (length(mR)>1)
    {
        for (i in 1:(length(mR)-1))
        {
            for (j in 2:length(mR))
            {
                Rij <- R[which[i]]*R[which[j]]
                if (Rij > 0)
                    start <- "+"
                else
                    start <- "-"
                Rij <- abs(Rij)
                Rij[Rij==1] <- ""
                sdevn2 <- paste(start, "2\\times ", Rij,
                                "Cov(\\hat{\\beta}_", which[i]-1,
                                ",\\hat{\\beta}_", which[j]-1,")",
                                sep="")
                sdevn <- paste(sdevn, sdevn2,sep="")
            }
        }                
    }
    what <- paste("Testing $H_0:~", coefnW, "=",q, "$ against $H_1:~", coefnW,
                  switch(alter, diff="\\neq", greater=">", less="<"),q, "$",
                  " at ", size*100, "\\%", sep="")    
    t <- m/s
    dist0 <- ifelse(exact, "\\sim ", "\\approx")
    dist <- ifelse(exact, paste("t_{",df,"}",sep=""), "N(0,1)")
    ftest <- paste("\\[\ntest = \\frac{",coefn,"}{se} = \\frac{",
                   round(m,digits),"}{",round(s,digits),"}=",
                   round(t,digits),dist0, dist,"\n\\]\n")
    ftest2 <- paste("where, \n\\[se=\\sqrt{",sdevn,"}\n\\]", sep="")
    ftest <- paste(ftest, ftest2)
    if (alter == "diff")
    {
        crit <- ifelse(exact, qt(1-size/2,df), qnorm(1-size/2))
        if (abs(t)>crit)
        {
            mes <- paste("Since $|", round(t,digits), "|>",
                         round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        } else {
            mes <- paste("Since $|", round(t,digits), "|<",
                         round(crit,digits),
                         "$ (the ", (1-size/2)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
        }
    } else if (alter == "greater") {
        crit <- ifelse(exact, qt(1-size,df), qnorm(1-size))
        if (t>crit)
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (1-size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    } else {
        crit <- ifelse(exact, qt(size,df), qnorm(size))
        if (t<crit)
            mes <- paste("Since $", round(t,digits), "<", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we reject $H_0$.", sep="")
        else
            mes <- paste("Since $", round(t,digits), ">", round(crit,digits),
                         "$ (the ", (size)*100, "\\% quantile of the $",
                         dist, "$), we do not reject $H_0$.", sep="")
    }
    if (!exact)
        mes <- paste(mes, " (The N(0,1) is an approximation based on the C.L.T. ",
                     "because either the distribution of the error term is unknown ",
                     "and/or the errors are not homoskedastic.)", sep="")    
    obj <- list(ftest=ftest, mes=mes,
                type="Testing linear combinations of regression coefficients",
                test=what)
    class(obj) <- "metricsSol"
    obj                            
}        

testRegF <- function(object, objectr=NULL, size=.05, digits=4,
                     assume = c("Normal", "nonNormal"), dround=NULL,
                     Rsq=NULL, SSRu=NULL, SSRr=NULL, nrest=NULL, n=NULL, ncoef=NULL)
{
    assume <- match.arg(assume)
    if ((!is.null(Rsq) | is.null(objectr)) & is.null(SSRu))
    {
        if (!is.null(Rsq))
        {
            if (any(c(is.null(n), is.null(ncoef))))
                stop("The sample size and/or the number of coefficients is missing")
            r2 <- Rsq
            df1 <- ncoef-1
            df2 <- n-ncoef
        } else {
            r2 <- summary(object)$r.squared
            if (!is.null(dround))
                r2 <- round(r2, dround)
            df1 <- length(coef(object))-1
            df2 <- object$df
        }
        if (assume=="Normal")
        {
            distAs <- "\\sim "
            dist <- paste("F(",df1,",",df2,")",sep="")
            div <- paste("/", df1, sep="")
            t <- r2/(1-r2)*df2/df1
        } else {
            distAs <- "\\approx "
            dist <- paste("\\chi^2_",df1,sep="")
            div <- ""
            t <- r2/(1-r2)*df2            
        }        
        ftest <- paste("\\[\ntest = \\frac{R^2",div,"}{(1-R^2)/",df2,
                       "} = \\frac{", round(r2,digits),div,
                       "}{(1-",round(r2,digits),")/",df2,"} = ",
                       round(t,digits), distAs, dist, "\n\\]\n")                
    } else {
        if (!is.null(SSRu))
        {
            if (any(c(is.null(SSRr),is.null(nrest), is.null(n), is.null(ncoef))))
                stop("For manual tests, the two SSR, the sample size, the number of constraints and the number of coefficients are needed")
            ssr_ur <- SSRu
            ssr_r <- SSRr
            df1 <- nrest
            df2 <- n-ncoef
        } else {
            ssr_ur <- sum(residuals(object)^2)
            ssr_r <- sum(residuals(objectr)^2)
            if (!is.null(dround))
            {
                ssr_ur <- round(ssr_ur,dround)
                ssr_r <- round(ssr_r,dround)
            }
            df1 <- objectr$df-object$df
            df2 <- object$df
        }
        if (assume=="Normal")
        {
            distAs <- "\\sim "
            dist <- paste("F(",df1,",",df2,")",sep="")
            div <- paste("/", df1, sep="")
            t <- (ssr_r-ssr_ur)/ssr_ur*df2/df1
        } else {
            distAs <- "\\approx "
            dist <- paste("\\chi^2_",df1,sep="")
            div <- ""
            t <- (ssr_r-ssr_ur)/ssr_ur*df2
        }        
        ftest <- paste("\\[\ntest = \\frac{(SSR_{r}-SSR_{ur})",
                       div,"}{SSR_{ur}/", df2,"}",
                       "= \\frac{(",round(ssr_r,digits),"-",
                       round(ssr_ur,digits),")",
                       div,"}{",round(ssr_ur,digits),"/", df2,"} = ",
                       round(t,digits), distAs, dist, "\n\\]\n")        
    }
    crit <- ifelse(assume=="Normal", qf(1-size,df1,df2), qchisq(1-size,df1))
    if (t>crit)
    {
        mes <- paste("Since $", round(t,digits), ">",
                     round(crit,digits),
                     "$ (the ", (1-size)*100, "\\% quantile of the $",
                     dist, "$), we reject $H_0$.", sep="")
    } else {
        mes <- paste("Since $", round(t,digits), "<",
                     round(crit,digits),
                     "$ (the ", (1-size)*100, "\\% quantile of the $",
                     dist, "$), we do not reject $H_0$.", sep="")
        
    }
    if (assume!="Normal")
        mes <- paste(mes, " (The $\\chi^2_{", df1,
                     "}$ is an approximation based on the C.L.T. ",
                     "because the distribution of the error term is unknown.", sep="")
    what <- "Testing joint hypotheses using the F-test"
    obj <- list(ftest=ftest, mes=mes, type=what, test=what)
    class(obj) <- "metricsSol"
    obj                            
}        


