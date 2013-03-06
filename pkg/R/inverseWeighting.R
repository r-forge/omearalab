require(pls)
#uses inverse weighting
#x0 is the vector containing the coordinates of a point to interpolate at (to do a grid, lapply this function)
#x is the sampled values (can be multivariate)
#y is the univariate response values (to do multivariate, call this multiple times with different y)
#centering is so that if the y values are on different scales, the ones smaller in magnitude aren't effectively ignored
#if you are wrapping this function, center y first and then do center=FALSE to speed up calculations
InverseInterpolateSingle <- function(x0, x, y, p, center=TRUE, scale=TRUE, dist.method="euclidean") {
	if(center+scale>0) {
		y<-stdize(y, center, scale)
	}
	distances<-(apply(y, 1, dist.mod, x=x, dist.method=dist.method))^p
	return(weighted.mean(y, w=distances)) 
}

dist.mod <- function(y, x, dist.method="euclidean") {
	return(dist(rbind(y,x), dist.method=dist.method))
}

FindP <- function(y, x, dist.method="euclidean") {
	#do cross validation to find estimate of p
	return(p)
}
