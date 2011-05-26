#takes a given consensus topology and assigns branch lengths based on the average of those source trees that have that branch
library(phylobase)

isClade<-function(phy,taxonlist) {
	if (length(taxonlist)==1) {
		return(TRUE)
	}
	if(length(descendants(phy,MRCA(phy,taxonlist),type="tips"))==length(taxonlist)) {
		return(TRUE)
	}
	else {
		return(FALSE)
	}
}

rootAndMakePhylo4<-function(phy) {
	return(as(root(phy,outgroup=c(1:length(phy$tip.label)),resolve.root=TRUE),"phylo4")
}

edgeLengthTaxset<-function(phy,taxonlist) {
	return(edgeLength(phy,MRCA(phy,taxonlist)))
}

summarizeNode<-function(nodeId,focalTree,sourceTreeList) {
	matchingVector<-unlist(lapply(sourceTreeList,isClade,descendants(focalTree,nodeId,type="tips")))
	proportion<-sum(matchingVector) / length(matchingVector)
	lengths<-unlist(lapply(sourceTreeList[matchingVector],edgeLengthTaxset,descendants(focalTree,nodeId,type="tips")))
	result<-c(nodeId,proportion,mean(lengths,na.rm=TRUE),median(lengths,na.rm=TRUE),sd(lengths,na.rm=TRUE))
	names(result)<-c("nodeId","proportion","mean_brlen","median_brlen","sd_brlen")
	print(result)
	return(result)
}

#trees should be rooted
consensusBrlen<-function(focalTree,sourceTreeList,type=c("proportion","mean_brlen","median_brlen","sd_brlen")) {
	type<-match.arg(type)
	if (class(focalTree)!="phylo4") {
		focalTree<-as(focalTree,"phylo4")
	}
	if (class(sourceTreeList[[1]])!="phylo4") {
		sourceTreeList<-lapply(sourceTreeList,as,"phylo4")
	}
	allNodes<-nodeId(focalTree,"all")
	allNodes<-allNodes[which(allNodes!=nodeId(focalTree,"root"))] #do not care about root edge
	allResults<-sapply(allNodes,summarizeNode,focalTree,sourceTreeList)
	newEdgeLengths<-edgeLength(focalTree)
	for (nodeIndex in 1:length(allNodes)) {
	#	print(paste("Old edge length for node",allNodes[nodeIndex],"is",edgeLength(focalTree,allNodes[nodeIndex])))
		print(names(newEdgeLengths))
		print(names(getEdge(focalTree,allNodes[nodeIndex]))[1])
		print(allResults[which(row.names(allResults)==type),nodeIndex])
		newEdgeLengths[ which(names(newEdgeLengths)==getEdge(focalTree,allNodes[nodeIndex])) ]<-allResults[which(row.names(allResults)==type),nodeIndex]
	#	print(paste("New edge length",edgeLength(focalTree,allNodes[nodeIndex])))
	#	print(focalTree)
	}
	edgeLength(focalTree)<-newEdgeLengths
#	print(allResults)
	return(focalTree)
}



##########Scratch space
sourceTreeList<-list(b)
for (i in 2:6) {
	sourceTreeList[[i]]<-b
}
for (i in 7:11) {
	sourceTreeList[[i]]<-as(root(rcoal(16),outgroup=c(1:16),resolve.root=TRUE),"phylo4")
}