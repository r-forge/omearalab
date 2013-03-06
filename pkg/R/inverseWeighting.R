require(nloptr)
#uses inverse weighting
#x.obs is the vector containing the coordinates of a point to interpolate at (to do a grid, lapply this function)
#x is the sampled values (can be multivariate)
#y is the univariate response values (to do multivariate, call this multiple times with different y)
#centering is so that if the y values are on different scales, the ones smaller in magnitude aren't effectively ignored
#if you are wrapping this function, center y first and then do center=FALSE to speed up calculations
InverseInterpolateSingle <- function(x.obs, x.mat, y.mat, p, standardize=TRUE, dist.method="euclidean") {
	if(standardize) {
	    x.mat.means <- colMeans(x.mat)
	    x.mat.sds <- apply(x.mat, 2, sd)
	    for (i in sequence(dim(x.mat)[2])) {
			x.mat[,i] <- (x.mat[,i] - x.mat.means[i]) / x.mat.sds[i]	    
			x.obs[i] <- (x.obs[i] - x.mat.means[i]) / x.mat.sds[i]
	    }
	}
	distances<-(apply(x.mat, 1, dist.mod, x.obs=x.obs, dist.method=dist.method))^p
	return(weighted.mean(y.mat, w=distances)) 
}

InverseInterpolateMany <- function(x.obsmat, x.mat, y.mat, p, standardize=TRUE, dist.method="euclidean") {
	if (standardize) {
	    x.mat.means <- colMeans(x.mat)
	    x.mat.sds <- apply(x.mat, 2, sd)
	    for (i in sequence(dim(x.mat)[2])) {
			x.mat[,i] <- (x.mat[,i] - x.mat.means[i]) / x.mat.sds[i]	    
			x.obsmat[,i] <- (x.obsmat[,i] - x.mat.means[i]) / x.mat.sds[i]
	    }
	}
	return(apply(x.obsmat, 1, InverseInterpolateSingle, x.mat=x.mat, y.mat=y.mat, p=p, standardize=FALSE, dist.method=dist.method))
}

dist.mod <- function(x.mat, x.obs, dist.method="euclidean") {
	return(dist(rbind(x.mat,x.obs), dist.method=dist.method))
}

FindP <- function(x.mat, y.mat, dist.method="euclidean") {
  opt <- nloptr(x0=1, eval_f=LOOV, lb=0.000001, ub=100, opts = list("algorithm"="NLOPT_LN_NEWUOA_BOUND"))
	return(opt$solution)
}

RMSE <- function(predictions, true) {
	return(sqrt(mean((predictions - true)^2)))
}

#p first to make optimization easier
LOOV <- function(p, x.mat, y.mat, dist.method="euclidean") {
	predictions <- sapply(sequence(dim(x.mat)[1]), x.mat=x.mat, y.mat=y.mat, p=p, dist.method=dist.method)
	return(RMSE(predictions, y.mat))
}

DropOne <- function(drop.index, x.mat, y.mat, p, dist.method="euclidean") {
	x.obs <- x.mat[drop.index, ]
	x.mat <- x.mat[-drop.index, ]
	return(InverseInterpolateSingle(x.obs, x.mat, y.mat, p, standardize=FALSE, dist.method=dist.method))
}