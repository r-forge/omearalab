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

colCountIf<-function(x,val) {
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

thetaindividual<-function(collapseMatrix,complete=FALSE,thetaMap) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,thetaMap=thetaMap)
	class(tI)<-"thetaindividual"
	return(tI)
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

kPopInterval<-function(popInterval) {
	#returns the number of free parameters needed for that interval object. For example, having everything collapse in one step requires one param (the tMRCA). Having one collapse then a second requires two. Having no collapse requires 0
	maxVector<-colMax(popInterval$collapseMatrix)
	return(length(which(maxVector>0)))
}

#the basic idea here is that at each population in each time interval there is a theta. These can all be set to the same value, allowed to vary, or assigned in clumps (i.e., pops 1, 3, and 6 have the same theta value)
#this generates all such mappings, subject to staying within the maximum number of free parameters
generateThetaIndividuals<-function(popVector,popIntervalsList=generateIntervals(popVector),maxK=max(1,floor(sum(popVector)/20))) {
	thetaIndividualsList<-list()
	for (i in 1:length(popIntervalsList)) {
		thetaMapTemplate<-1+0*popIntervalsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		numLineages=sum(thetaMapTemplate,na.rm=TRUE)
		possibleMappings<-compositions(numLineages)
		for (mappingIndex in 1:dim(possibleMappings)[2]) {
			thisMapping<-possibleMappings[,mappingIndex]
			if ((length(which(thisMapping>0))+kPopInterval(popIntervalsList[[i]]) )<=maxK) { #only do it on those mappings that have not too many free parameters
				thetaMap<-thetaMapTemplate
				whichPositions <- which(thetaMap==1)
				for (positionIndex in 1:length(whichPositions)) {
					position=whichPositions[positionIndex]
					paramPosition<-which(thisMapping>0)[1]
					thetaMap[position]=paramPosition #the position of the first parameter
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we've used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				thetaIndividualsList[[length(thetaIndividualsList)+1]]<-thetaindividual(popIntervalsList[[i]]$collapseMatrix, popIntervalsList[[i]]$complete, thetaMap)
			}
		}
	}
	return(thetaIndividualsList)
}

#now we will generate all possible assignments of pairwise migration. Again, we want to keep the total number of free parameters (times, thetas, migration rates) under our chosen max
generateMigrationIndividuals<-function(popVector,thetaIndividualsList=generateThetaIndividuals(popVector), maxK=max(1,floor(sum(popVector)/20))) {
	migrationIndividualsList<-list()
	for (i in 1:length(thetaIndividualsList)) {
		#things to consider: migration between each subpop, as subpops coalesce
	}
}
