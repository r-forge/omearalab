library(partitions)

colMax<-function(x,na.rm=TRUE) {
	maxVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		maxVal[i]<-max(x[,i],na.rm=na.rm)
	}
	return(maxVal)
}

colMin<-function(x,na.rm=TRUE) {
	minVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		minVal[i]<-min(x[,i],na.rm=na.rm)
	}
	return(minVal)
}

generateIntervals<-function(nPop,partition,intervalList) {
	
}


generateScenarios<-function(popVector,maxK=max(1,floor(sum(popVector)/20)) {
	#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
	#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
	nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,colMax(firstIntervals)==1]
	firstIntervals<-firstIntervals[,colSums(firstIntervals)>1] #which populations join in the first interval (counting back from the tips)
	#TODO: Now have to walk back recursively, getting all intervals possible given each of the initial intervals
	#probably use recursive fn
}


library(phylobase)
library(phangorn)
#old code
generateScenarios<-function(popVector,maxK=max(1,floor(sum(popVector)/20)),maxMigration=2,maxSplit=2,maxBottleneck=2) {
	#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
	#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
	nPop <- length(popVector)
	if (nPop<2) {
		stop("assumes at least two populations, you have ",nPop)
	}
	allHistories<-allTrees(nPop,rooted=TRUE,tip.label=c(1:nPop)) #all bifurcating histories


	#TODO: now to calculate all multifurcating trees, too. These are also stored in allHistories
	
	for (treeIndex in 1:length(allHistories)) {
		focalTree<-as(allHistories[[treeIndex]],"phylo4")
		
		#now have to do all possible node orderings: for ((a,b),(c,d)), which pair coalesces first? if a and b, then can have migration from ab_ancestor to c or to d, but not if reversed
		orderedTrees<-list() #TODO
		
		for (otreeIndex in 1:length(orderedTrees)) {
			orderedTree<-orderedTrees[[otreeIndex]]
			
		}
	}
}