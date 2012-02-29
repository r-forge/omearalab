#a function to deal with Gould-style punctuated equilibrium: at a speciation event ONE of the two descendant lineage changes. Figuring out which one looks NP-nasy. So be it.

#do this with phylo, rather than phylo4, as we'll be playing with this with geiger fns for some

library(geiger)
library(sfsmisc)
library(partitions) #for converting from binary back to decimal
library(gmp) #for dealing with big integers
library(optimx)
library(phylobase)
library(rgenoud)
library(snow)


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
#	print(paste("transformation.param=",transformation.param))
	if(data.type=="Continuous") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitContinuous(phy,data)[[1]]$lnl) #want to minimize neg lnL
		#print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newLnL
		}
	}
	if(data.type=="Discrete") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitDiscrete(phy,data,model=data.model)[[1]]$lnl)
		#print(newNegLnL)
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
 #print(transformation.param)
	
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
 #print(summary(phyClock))
 #print(summary(phyGould))
	
	phyClock<-transformation.fn(phyClock,stretchParam)
  #print(summary(phyClock))

	phyClockHeight<-max(depthTips(as(phyClock,"phylo4")))
	phyGouldHeight<-max(depthTips(as(phyGould,"phylo4")))
#	if(sum(sum(phyClock$edges[,1]==phyGould$edges[,1]),sum(phyClock$edges[,2]==phyGould$edges[,2]))!=2*max(dim(phyClock$edges)[1],dim(phyGould$edges)[1])) {
#    print(cbind(phyGould$edges,phyClock$edges))
#    stop("Imperfect match between phyClock and phyGould: they must match exactly, including internal ordering of information")
#	}
	
	phy<-phyClock
	phy$edge.length<-(gouldWeight) * (phyGould$edge.length/phyGouldHeight) + (1 - gouldWeight) * (phyClock$edge.length/phyClockHeight)
  plot(phy)
 neglnL<-badVal
	if(data.type=="Continuous") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitContinuous(phy,data)[[1]]$lnl) #want to minimize neg lnL
		#print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newNegLnL
		}
	}
	if(data.type=="Discrete") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitDiscrete(phy,data,model=data.model)[[1]]$lnl)
		#print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newNegLnL
		}
	}
	return(neglnL)
}

likelihoodGouldPlusNonGouldTransformPlusError<-function(transformation.param,phyClock,phyGould,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,ln.transformed=FALSE) {
  data.type<-match.arg(data.type)
  if(ln.transformed) {
    transformation.param<-exp(transformation.param) 
  }
	gouldWeight<-transformation.param[1]
	stretchParam<-transformation.param[2]
  errorParam<-transformation.param[3]
  print(transformation.param)
	print("Starting bounds checking")
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
  if (errorParam<0) {
    return(badVal) 
  }
  print("Finished bounds checking")
 #print(summary(phyClock))
 #print(summary(phyGould))
	
	phyClock<-transformation.fn(phyClock,stretchParam)
  #print(summary(phyClock))

	phyClockHeight<-max(depthTips(as(phyClock,"phylo4")))
	phyGouldHeight<-max(depthTips(as(phyGould,"phylo4")))
#	if(sum(sum(phyClock$edges[,1]==phyGould$edges[,1]),sum(phyClock$edges[,2]==phyGould$edges[,2]))!=2*max(dim(phyClock$edges)[1],dim(phyGould$edges)[1])) {
#    print(cbind(phyGould$edges,phyClock$edges))
#    stop("Imperfect match between phyClock and phyGould: they must match exactly, including internal ordering of information")
#	}
	
	phy<-phyClock
	phy$edge.length<-(gouldWeight) * (phyGould$edge.length/phyGouldHeight) + (1 - gouldWeight) * (phyClock$edge.length/phyClockHeight)
  phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]+errorParam
  #plot(phy)
 neglnL<-badVal
  sink("/dev/null")
	if(data.type=="Continuous") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitContinuous(phy,data)[[1]]$lnl) #want to minimize neg lnL
		#print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newNegLnL
		}
	}
	if(data.type=="Discrete") {
		newNegLnL<-NA
		try(newNegLnL<-(-1)*fitDiscrete(phy,data,model=data.model)[[1]]$lnl)
		#print(newNegLnL)
		if(is.finite(newNegLnL)) {
			neglnL<-newNegLnL
		}
	}
  sink()
  print(paste("negLnL= = ",neglnL))

	return(neglnL)
}


fitNonGouldTransform<-function(phy,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,optimx.method="Nelder-Mead") {
	data.type<-match.arg(data.type)
	results<-optimx(par=c(1),fn=likelihoodNonGouldTransform,method=optimx.method,phy=phy,data=data,transformation.fn=transformation.fn,data.type=data.type,data.model=data.model,badVal=badVal)
	return(results)
}

fitGouldPlusNonGouldTransform<-function(phy,data,gouldVector,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,optimx.method="Nelder-Mead") {
  data.type<-match.arg(data.type)
  results<-optimx(par=c(0.5,1),fn=likelihoodGouldPlusNonGouldTransform,method=optimx.method,phyClock=phy,phyGould=transformGould(phy,gouldVector),data=data,transformation.fn=transformation.fn,data.type=data.type,data.model=data.model,badVal=badVal)
	return(results)
}

fitGouldPlusNonGouldPlusErrorTransform<-function(phy,data,gouldVector,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,optimx.method="Nelder-Mead") {
  data.type<-match.arg(data.type)
  results<-optimx(par=c(0.5,1,0.0000000001),fn=likelihoodGouldPlusNonGouldTransformPlusError,method=optimx.method,phyClock=phy,phyGould=transformGould(phy,gouldVector),data=data,transformation.fn=transformation.fn,data.type=data.type,data.model=data.model,badVal=badVal)
  return(results)
}

likGouldPlusNonGouldPlusErrorTransformContinuousParamOptimized<-function(gouldVector,phy,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,optimx.method="Nelder-Mead",itnmax=NULL) {
  data.type<-match.arg(data.type)
  print(gouldVector)
  if(min(gouldVector)<0) {
    print("Gould vector min too low")
    return(badVal) 
  }
  if(max(gouldVector)>1) {
    print("Gould vector max too high")
    return(badVal) 
  }
  results<-optimx(par=log(c(0.5,1,0.0000000001)),fn=likelihoodGouldPlusNonGouldTransformPlusError,method=optimx.method,phyClock=phy,phyGould=transformGould(phy,gouldVector),data=data,transformation.fn=transformation.fn,data.type=data.type,data.model=data.model,badVal=badVal,ln.transformed=TRUE,itnmax=itnmax)
  print(results$fvalues)
  return(results$fvalues)
}

fitLikGouldPlusNonGouldPlusErrorTransformContinuousParamOptimized<-function(phy,data,transformation.fn,data.type=c("Continuous","Discrete"),data.model="ER",badVal=1000000000,optimx.method="Nelder-Mead",pop.size=1000,itnmax=NULL) {
  library(optimx)
  data.type<-match.arg(data.type)
  starting.values<-matrix(data=rbinom(pop.size*lengthGouldVector(phy),1,0.5),nrow=pop.size)
  Domains<-matrix(c(rep(0,lengthGouldVector(phy)),rep(1,lengthGouldVector(phy))),ncol=2,byrow=FALSE)
  results<-genoud(fn=likGouldPlusNonGouldPlusErrorTransformContinuousParamOptimized,nvars=lengthGouldVector(phy),boundary.enforcement=2, max=FALSE,data.type.int=TRUE,pop.size=pop.size,starting.values=starting.values, Domains=Domains,            phy=phy,data=data,transformation.fn=transformation.fn,data.type=data.type,data.model=data.model,badVal=badVal,itnmax=itnmax)
  return(results)
}


#Add measurement error: important with zero length branches

################## test data ###################

phy<-rcoal(8)
trueGouldVector<-toGouldVector(5,2,lengthGouldVector(phy))
phyGouldTrue<-nonzeroUniformTree(transformGould(phy,trueGouldVector))
data<-sim.char(phyGouldTrue,matrix(1),1)[,,1]
actualGouldVector<-trueGouldVector

#result<-likelihoodGouldPlusNonGouldTransform(c(1,1),phy,nonzeroUniformTree(transformGould(phy,actualGouldVector)),data,kappaTree,"Continuous")

full<-fitLikGouldPlusNonGouldPlusErrorTransformContinuousParamOptimized(phy,data,kappaTree,"Continuous",pop.size=100,itnmax=30)


liks<-c()
liks.1<-c()
distances<-c()
gouldVectors<-matrix(trueGouldVector,nrow=1)
for (i in as.numeric(minIndexGouldVector(phy)):as.numeric(maxIndexGouldVector(phy))) {
  plot(transformGould(phy,toGouldVector(i,2,lengthGouldVector(phy))))
  phyClock<-phy
  distances<-append(distances,sum(abs(actualGouldVector-toGouldVector(i,2,lengthGouldVector(phy)))))
  phyGould<-nonzeroUniformTree(transformGould(phy,toGouldVector(i,2,lengthGouldVector(phy))))
  #phyClock<-kappaTree(phyClock,runif(1))
  liks<-append(liks,likelihoodGouldPlusNonGouldTransformPlusError(c(1,1,0.000001),phyClock,phyGould,data,kappaTree,"Continuous"))
  liks.1<-append(liks.1,likelihoodGouldPlusNonGouldTransformPlusError(c(1,1,.1),phyClock,phyGould,data,kappaTree,"Continuous"))
  #Sys.sleep(5)
  gouldVectors<-rbind(gouldVectors,toGouldVector(i,2,lengthGouldVector(phy)))
}
gouldVectors<-gouldVectors[-1,]
plot(liks,liks.1,type="n")
text(liks,liks.1,distances)
regression<-lm(liks~gouldVectors)
