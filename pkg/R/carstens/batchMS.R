library(partitions)

colMax<-function(x,na.rm=TRUE) {
	print("colMax")
	maxVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		maxVal[i]<-max(x[,i],na.rm=na.rm)
	}
	return(maxVal)
}

colMin<-function(x,na.rm=TRUE) {
	print("colMin")
	minVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		minVal[i]<-min(x[,i],na.rm=na.rm)
	}
	return(minVal)
}

colCountIf<-function(x,val) {
	print("colCountIf")
	print(paste("dim(x)=",dim(x)))
	countVec<-rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		countVec[i]=length(which(x[,i]==val))
	}
	return(countVec)
}

colCountLast<-function(x,val) {
	return(length(which(x[,dim(x)[2]]==val)))
}



popinterval<-function(collapseMatrix,complete=FALSE) {
	#collapse matrix has populations as rows and each generation going back in time as columns
	pInt<-list(collapseMatrix=collapseMatrix,complete=complete)
	class(pInt)="popinterval"
	if (max(collapseMatrix,na.rm=TRUE)==0) {
		pInt$complete=TRUE #if the model is one of no collapse, nothing further to do
	}
	return(pInt)
}

returnIncompletes<-function(popIntervalsList) {
	incompleteElements<-c()
	for (i in 1:length(popIntervalsList)) {
		if (!popIntervalsList[[i]]$complete) {
			incompleteElements<-append(incompleteElements,i)
		}
	}
	return(incompleteElements)
}

updateCompletes<-function(popIntervalsList) {
	for (i in 1:length(popIntervalsList)) {
		if (colCountLast(popIntervalsList[[i]]$collapseMatrix,1)==0) {
			popIntervalsList[[i]]$complete<-TRUE #end up with having no more collapses
		}
		else if (colCountLast(popIntervalsList[[i]]$collapseMatrix,0)==0) {
			popIntervalsList[[i]]$complete<-TRUE #end up as one population
		}
	}
	return(popIntervalsList)	
}

completeIntervals<-function(popIntervalsList) {
	while(length(returnIncompletes(popIntervalsList))>0) {
		incompleteElements=returnIncompletes(popIntervalsList)
		newPopIntervalsList<-popIntervalsList[-1*incompleteElements] #has only complete elements
		for (index in 1:length(incompleteElements)) {
			currentParent=popIntervalsList[[incompleteElements[index] ]]
			lastGen<-currentParent$collapseMatrix[,dim(currentParent$collapseMatrix)[2] ]
			#all those in the lastGen in state 1 were merged into one population, id'ed by the one with the lowest id
			survPop=length(which(lastGen==0))
			if (length(which(lastGen==1))>0) {
				survPop=survPop+1
			}
			else {
				stop("this population was complete")
			}
			rawIntervals<-c()
			if (survPop>1) {
				rawIntervals<-blockparts(c(1:survPop),survPop,include.fewer=TRUE)
				rawIntervals<-rawIntervals[,colMax(rawIntervals)<=1] #we're okay with having all zeros: no population collapse
				rawIntervals<-rawIntervals[,colCountIf(rawIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
			}
			else {
				rawIntervals=matrix(c(0,1),nrow=1)
			}
			rowMapping<-c(min(which(lastGen==1)),which(lastGen==0))
			rawIntervalsRescaled<-matrix(NA,ncol=dim(rawIntervals)[2],nrow=length(lastGen))
			for (j in 1:length(rowMapping)) {
				rawIntervalsRescaled[rowMapping[j], ]<-rawIntervals[j,]
			}
			for (k in 1:dim(rawIntervalsRescaled)[2]) {
				newPopIntervalsList[[length(newPopIntervalsList)+1]]<-popinterval(cbind(currentParent$collapseMatrix,matrix(rawIntervalsRescaled[,k],ncol=1)))
			}
		}
		popIntervalsList<-newPopIntervalsList
		popIntervalsList<-updateCompletes(popIntervalsList)
	}
	return(popIntervalsList)
}


generateIntervals<-function(popVector,maxK=max(1,floor(sum(popVector)/20))) {
	#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
	#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
	nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,colMax(firstIntervals)<=1] #we're okay with having all zeros: no population collapse
	firstIntervals<-firstIntervals[,colCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	popIntervalsList<-list()
	for (i in 1:dim(firstIntervals)[2]) {
		popIntervalsList[[i]]<-popinterval(as.matrix(firstIntervals[,i],ncol=1))
	}
	popIntervalsList<-completeIntervals(updateCompletes(popIntervalsList))
	return(popIntervalsList)
}

