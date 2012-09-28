
library(laser)

##Input is phy, r=net diversification, e=extinction fraction, k=carrying capacity
test.laser.lik<-function(phy, r, e ,k){
	
	x <- branching.times(phy)
	
	obj<-NULL
	obj$bd.LH <- calcLHbd(x, r, e)
	obj$DD.bd.LH <- DDL(x, r, k)$LH
	obj$diff.DD.bd_bd <- obj$DD.bd.LH - obj$bd.LH
	
	obj
}


DDL <- function(x, r, k){
	if (!is.numeric(x)) stop("object x not of class 'numeric'")
	#calculates likelihoods under DD logistic model
	x <- rev(sort(x))
	N <- length(x)+1
	b <- sort(x)
	z <- rev(c(b[1], diff(b)))
	r <- r
	k <- k
	
	ddfunc <- function(r, k){
		-(sum(log(2:(N-1))) + (N-2)*log(r) + sum(log(1-((2:(N-1))/k))) - sum((2:N)*r*z) + sum(z*r*(2:N)^2)/k)
	}
	res <- ddfunc(r, k)
	#may want to recode this to use 'optim' rather than 'nlm' -- Dan, you do not (JMB)
	aic <- 2*res + 4
	summ <- structure(list(LH = -res, aic = aic, r1 = r, kparam = k))
	return(summ)
}

