Quadra <- function(a,b,c)
	{
	if(a==0)
		stop("It is not a quadratique function;
				 'a' must be different from zero")
	obj <- list(a=a,b=b,c=c)
	class(obj) <- "Quadra"
	return(obj)
	}

print.Quadra <- function(x, ...)
	{
	cat("\nSecond order polynomial\n\n")
	cat("F(x) = Ax^2 + Bx + C\n")
	cat("with: A=",x$a,", B=", x$b, ", C=", x$c,"\n\n")	
	}

zeros <- function(object, ...)
  {
  UseMethod("zeros")
  }

zeros.Quadra <- function(object, ...)
	{
	det <- object$b^2-4*object$a*object$c
	if (det>.Machine$double.eps)
		{
		r1 <- (-object$b-sqrt(det))/(2*object$a)
		r2 <- (-object$b+sqrt(det))/(2*object$a)
		r <- cbind(r1,r2)
		class(r) <- "zeros"
		attr(r,"type") = "Real and distinct"
		}
	if (abs(det) <= .Machine$double.eps)
		{
		r1 <- -object$b/(2*object$a)
		r <- cbind(r1,r1)
		class(r) <- "zeros"
		attr(r,"type") = "Real and identical"
		}
	if (det < -.Machine$double.eps)
		{
		det <- sqrt(-det)/(2*object$a)
		r1 <- -object$b/(2*object$a) - det*1i
		r2 <- -object$b/(2*object$a) + det*1i
		r <- cbind(r1,r2)
		class(r) <- "zeros"
		attr(r,"type") = "Complexe"
		}
	return(r)
	}

print.zeros <- function(x, ...)
	{
	n <- length(x)
	cat("\nType of zeros: ", attr(x,"type"),"\n\n")
	for (i in 1:n)
		cat("Zero[",i,"] = ",x[i],"\n")
	cat("\n")
	}


solveP <- function(obj, ...)
  {
  UseMethod("solveP")
  }

solveP.Quadra <- function(obj, ...)
		{
		x <- -obj$b/(2*obj$a)
		f <- obj$a*x^2+obj$b*x+obj$c
		if (obj$a>0)
			what <- "min"
		else
			what <- "max"
		ans <- list(x=x,f=f,what=what)
		class(ans) <- "solveP.Quadra"
		return(ans)
		}

print.solveP.Quadra <- function(x, ...)
	{
	if (x$what=="min")
		mes <- "\nThe polynomial has a minimum at "
	else
		mes <- "\nThe polynomial has a maximum at "
	cat(mes,"x = ", x$x,"\n")
	cat("At that point, f(x) = ",x$f,"\n\n")
	}

addQuadra <- function(Q1,Q2)
		{
		if (class(Q1)!="Quadra" | class(Q2)!="Quadra")
			stop("This operator can only be applied to
				objects of class Quadra")	
		a <- Q1$a+Q2$a
		b <- Q1$b+Q2$b
		c <- Q1$c+Q2$c
		Quadra(a,b,c)
		}

"%+%" <- function(Q1,Q2)
	addQuadra(Q1,Q2)

plot.Quadra <- function(x,from=NULL,to=NULL, ...)
 	{
	f <- function(y)
		x$a*y^2+x$b*y+x$c

	res <- solveP(x)

	if(is.null(from) | is.null(to))
 		{
		from <- res$x-4
		to <- res$x+4
		}
	if (res$what=="min")
		{
		d <- max(f(to),f(from))-res$f 
		mes <- paste("Min=(",round(res$x,2),", ",round(res$f,2),")",sep="")
		}
	if (res$what=="max")
		{
		mes <- paste("Max=(",round(res$x,2),", ",round(res$f,2),")",sep="")
		d <- res$f - min(f(to),f(from))
		}
	
	curve(f,from,to,xlab="X",ylab="f(X)")
	if (x$b>0 & x$c>0)
		title(substitute(f(X)==a*X^2+b*X+c,x))
	if (x$b<0 & x$c>0)
		title(substitute(f(X)==a*X^2-b2*X+c,c(x,b2=-x$b)))
	if (x$b>0 & x$c<0)
		title(substitute(f(X)==a*X^2+b*X-c2,c(x,c2=-x$c)))
	if (x$b==0 & x$c>0)
		title(substitute(f(X)==a*X^2+c,x))
	if (x$b==0 & x$c<0)
		title(substitute(f(X)==a*X^2-c2,c(x,c2=-x$c)))
	if (x$c==0 & x$b>0)
		title(substitute(f(X)==a*X^2+b*x,x))
	if (x$c==0 & x$b<0)
		title(substitute(f(X)==a*X^2-b2*x,c(x,b2=-x$b)))

	points(res$x,res$f,col=3,cex=.8, pch=21,bg=3)
	if(res$what=="min")
		{
		text(res$x,res$f+.2*d,mes)
		arrows(res$x,res$f+.18*d,res$x,res$f)
		}
	else
		{
		text(res$x,res$f-.2*d,mes)
		arrows(res$x,res$f-.18*d,res$x,res$f)
		}

	z <- zeros(x)
	if (attr(z,"type")=="Real and distinct")
		{
		points(z[1],0,col=2,cex=.8,pch=21,bg=2)
		points(z[2],0,col=2,cex=.8,pch=21,bg=2)
		r1 <- paste(round(min(z),2))
		r2 <- paste(round(max(z),2))
		if(res$what=="min")
			{
			if(abs(res$f)>d/2)
				d2 <- -d
			else
				d2 <- d
			text(min(z),.25*d2,r1)
			text(max(z),.25*d2,r2)
			arrows(min(z),.23*d2,min(z),0)
			arrows(max(z),.23*d2,max(z),0)	
			}
		else
			{
			if(abs(res$f)>d/2)
				d2 <- -d
			else
				d2 <- d
			text(min(z),-.25*d2,r1)
			text(max(z),-.25*d2,r2)
			arrows(min(z),-.23*d2,min(z),0)
			arrows(max(z),-.23*d2,max(z),0)	
			}
		}
	if(attr(z,"type")!="Complexe" | attr(z,"type") == "Real and identical")
		abline(h=0)
	 }

summary.Quadra <- function(object, ...)
	{
	print(object)
	print(zeros(object))
	print(solveP(object))
	}
