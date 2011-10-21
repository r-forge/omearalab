#a function to deal with Gould-style punctuated equilibrium: at a speciation event ONE of the two descendant lineage changes. Figuring out which one looks NP-nasy. So be it.

#do this with phylo, rather than phylo4, as we'll be playing with this with geiger fns for some

library(geiger)
library(sfsmisc)
library(partitions) #for converting from binary back to decimal
library(gmp) #for dealing with big integers
library(optimx)
library(phylobase)


lengthGouldVector<-function(phy) {
	if(is.binary.tree(phy)) {
		return(Nnode(phy,internal.only=TRUE))
	}
	else {
		warning("Tree for Gould transform must be binary")
		return(NA)
	}
}

#utility functions
toGouldVector<-function (x, base = 2, ndigits = 1 + floor(log(max(x),base)))  #modification from sfsmisc package's digitsBase to deal with large numbers
{
	if (class(x)!="bigz") {
		x<-as.bigz(x)
	}
	if (class(base)!="bigz") {
		base<-as.bigz(base)
	}
	binaryVector<-c()
	while (x>0) {
		binaryVector<-c(as.numeric(x%%base),binaryVector)
		x<-x%/%base
	}
	fullVector<-rep(0,ndigits)
	lengthDiff<-length(fullVector)-length(binaryVector)
	if (lengthDiff<0) {
		lengthDiff<-0
		fullVector<-rep(0,length(binaryVector))
	}
	if (length(binaryVector)>0) {
		for (i in 1:length(binaryVector)) {
			fullVector[i+lengthDiff]<-binaryVector[i]
		}
	}
	return(fullVector)
}

minIndexGouldVector<-function(phy) {
	return(as.bigz(0)) #yeah, I know. But makes interface easier than having to remember this
}

maxIndexGouldVector<-function(phy) {
	return(-1+pow.bigz(2,lengthGouldVector(phy))) #why bigz? because this can get huge
}

transformGould<-function(phy,gouldVector) {
	brlenMatrix<-cbind(phy$edge,phy$edge.length)
	minInternalNode<-min(brlenMatrix[,1])
	maxInternalNode<-max(brlenMatrix[,1])
	for (focalNode in minInternalNode:maxInternalNode) {
		focalRowSet<-which(brlenMatrix[,1]==focalNode)
		whichElement<-1+gouldVector[focalNode-minInternalNode+1]
		focalRow<-focalRowSet[whichElement]
		brlenMatrix[focalRow,3]<-0
	}
	phy$edge.length<-brlenMatrix[,3]
	return(phy)
}

noTransform<-function(phy,...) {
	return(phy)
}

nonzeroUniformTree<-function(phy,...) {
	phy$edge.length[which(phy$edge.length>0)]<-1
	return(phy)
}

#this takes a tree, some data, and a transformation function that uses a single parameter (kappaTree or something similar). 
#it returns the optimal values of the parameters and the likelihood
likelihoodNonGouldTransform<-function(transformation.param,phy,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000) {
	data.type<-match.arg(data.type)
	phy<-transformation.fn(phy,transformation.param)
	neglnL<-badVal
	print(paste("transformation.param=",transformation.param))
	if(data.type=="Continuous") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitContinuous(phy,data)[[1]]$lnl) #want to minimize neg lnL
		print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newLnL
		}
	}
	if(data.type=="Discrete") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitDiscrete(phy,data,model=data.model)[[1]]$lnl)
		print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newNegLnL
		}
	}
	return(neglnL)
}

likelihoodGouldPlusNonGouldTransform<-function(transformation.param,phyClock,phyGould,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000) {
	data.type<-match.arg(data.type)
	gouldWeight<-transformation.param[1]
	stretchParam<-transformation.param[2]
	
	#keep in bounds
	if (gouldWeight<0) {
		return(badVal)
	}
	if (gouldWeight>1) {
		return(badVal)
	}
	if (stretchParam<0) {
		return(badVal)
	}
	if (stretchParam>1) {
		return(badVal)
	}
	
	phyClock<-transformation.fn(phyClock,stretchParam)
	phyClockHeight<-max(depthTips(as(phyClock,"phylo4")))
	phyGouldHeight<-max(depthTips(as(phyGould,"phylo4")))
	
	if(sum(sum(phyClock$edges[,1]==phyGould$edges[,1]),sum(phyClock$edges[,2]==phyGould$edges[,2]))!=2*max(dim(phyClock$edges)[1],dim(phyGould$edges)[1])) {
		stop("Imperfect match between phyClock and phyGould: they must match exactly, including internal ordering of information")
	}
	
	phy<-phyClock
	phy$edge.length<-(gouldWeight) * (phyClock$edge.length/phyGouldHeight) + (1 - gouldWeight) * (phyClock$edge.length/phyClockHeight)
	neglnL<-badVal
	if(data.type=="Continuous") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitContinuous(phy,data)[[1]]$lnl) #want to minimize neg lnL
		print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newLnL
		}
	}
	if(data.type=="Discrete") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitDiscrete(phy,data,model=data.model)[[1]]$lnl)
		print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newNegLnL
		}
	}
	return(neglnL)
}

fitNonGouldTransform<-function(phy,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,optimx.method="Nelder-Mead") {
	data.type<-match.arg(data.type)
	results<-optimx(par=c(1),fn=likelihoodNonGouldTransform,method=optimx.method,phy=phy,data=data,transformation.fn=transformation.fn,data.type=data.type,data.model=data.model,badVal=badVal)
	return(results)
}



