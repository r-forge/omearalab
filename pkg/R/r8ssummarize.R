#this takes a single tree that has been run through r8s to get node ages and also a set of bootstrapped trees (same topology, different characters and input brlen).
#it then plots the single tree and includes node labels that are the 95% CI from the bootstrap trees
#one nice feature is that it works even if the input trees differ in which nodes have been collapsed (in fact, it should work with completely different bootstrap tree topologies, but haven't tried that)
edr8sSummarize<-function() {
	setwd("/Users/bomeara/Dropbox/CollabSchilling/consensusRun")
	
	library(ape)
	library(phylobase)
	besttree<-as(read.nexus(file="mp.con.brlen.r8s.FinalTrees"),"phylo4")
	boottrees<-read.nexus(file="treeboot.brlen.tre.r8s.FinalTrees")
	#boottrees<-read.nexus(file="mp.con.brlen.r8s.FinalTrees")
	newNodeLabels<-nodeLabels(besttree)
	nodeNames<-names(newNodeLabels)
	print(length(boottrees))
	multipliers=c(1,92.4/77)
	overallheights=c(77,92.4)
	
	for (loop in 1:2) {
	multiplier=multipliers[loop]
	overallHeight=overallheights[loop]
	for (nodeIndex in nNodes(besttree):1) {
			nodeDescendants=descendants(besttree,as.numeric(nodeNames[nodeIndex]),type="tips")
			#print(paste("node: ",as.numeric(nodeNames[nodeIndex])," descendants tips: "))
			#print(nodeDescendants)
			#print(" descendants all: ")
			#print(descendants(besttree,as.numeric(nodeNames[nodeIndex]),type="all"))
			#print(" ")
			heights=rep(NA,length(boottrees))
			print(nodeIndex)
	
			for (bootIndex in 1:length(boottrees)) {
				boottree<-as(boottrees[[bootIndex]],"phylo4")
				mrcaNode<-as.numeric(MRCA(boottree,nodeDescendants))
				
					#print(mrcaNode)
					#print(rootNode(boottree))
					#print(as.numeric(sumEdgeLength(boottree,mrcaNode)))
					heights[bootIndex]<-overallHeight-(multiplier*sumEdgeLength(boottree,ancestors(boottree,mrcaNode,type="ALL")))
				
			}
			#print(heights)
			if (as.numeric(nodeNames[nodeIndex])!=rootNode(besttree)) {
				newNodeLabels[nodeIndex]<-paste(sprintf("%.1f",overallHeight-(multiplier*sumEdgeLength(besttree,ancestors(besttree,as.numeric(nodeNames[nodeIndex]),type="ALL"))))," [",sprintf("%.1f",quantile(heights,0.025)),"-",sprintf("%.1f",quantile(heights,0.975)),"]",sep="")
			}
	}
	nodeLabels(besttree)<-newNodeLabels
	quartz()
	besttree.phylo<-as(besttree,"phylo")
	besttree.phylo$edge.length<-besttree.phylo$edge.length*multiplier
	plot.phylo(besttree.phylo)
	nodelabels(besttree.phylo$node.label,adj=c(0,0.5),bg="white",frame="none",cex=0.6)
	axisPhylo()
	}
}