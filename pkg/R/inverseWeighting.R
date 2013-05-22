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
	return(weighted.mean(y.mat, w=1/distances)) 
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

dist.mod <- function(x.mat, x.obs, dist.method="euclidean", standardize=FALSE) {
	if(standardize) {
	    x.mat.means <- colMeans(x.mat)
	    x.mat.sds <- apply(x.mat, 2, sd)
	    for (i in sequence(dim(x.mat)[2])) {
			x.mat[,i] <- (x.mat[,i] - x.mat.means[i]) / x.mat.sds[i]	    
			x.obs[i] <- (x.obs[i] - x.mat.means[i]) / x.mat.sds[i]
	    }
	}
	return(dist(rbind(x.mat,x.obs), method=dist.method))
}

FindP <- function(x.mat, y.mat, dist.method="euclidean") {
  return(optimize(f=LOOV, interval=c(0.00001, 100), x.mat=x.mat, y.mat=y.mat, dist.method="euclidean")$minimum)
}

RMSE <- function(predictions, true) {
	return(sqrt(mean((predictions - true)^2)))
}

#p first to make optimization easier
LOOV <- function(p, x.mat, y.mat, dist.method="euclidean") {
	predictions <- sapply(sequence(dim(x.mat)[1]), DropOne, x.mat=x.mat, y.mat=y.mat, p=p, dist.method=dist.method)
	return(RMSE(predictions, y.mat)) #jeremy made me do it this way
}

DropOne <- function(drop.index, x.mat, y.mat, p, dist.method="euclidean") {
  ncolx <- dim(x.mat)[2]
	x.obs <- matrix(x.mat[drop.index, ], ncol=ncolx)
	x.mat <- matrix(x.mat[-drop.index, ], ncol=ncolx)
  y.mat <- matrix(y.mat[-drop.index, ], ncol=dim(y.mat)[2])
	return(InverseInterpolateSingle(x.obs=x.obs, x.mat=x.mat, y.mat=y.mat, p=p, standardize=FALSE, dist.method=dist.method))
}