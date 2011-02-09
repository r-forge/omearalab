#this is a function to convert from phylo or phylo4 to a newick string, by default in simmap1.1 format. 
#By giving it a vector of taxon labels, it can give everything in the smallest clade containing those taxa
#a simmap state of 1, everything else a state of 0, and the edge leading to that clade can be broken up: 
#by default, half is in each state, but if change.position=0.75, the first 3/4 is in state 0 and the second 1/4
#is in state 1.
#If that edge is set to be entirely in one state (change.position==1 or change.position==0), there will be 
#one state with zero length on that branch. This can be optionally deleted using suppress.zero=TRUE

library(phylobase)
generate.simmap<-function(x, taxa.vector, change.position=0.5, digits=10, suppress.zero=FALSE,format="simmap1.1") {
	if(class(x)!="phylo4") {
		x<-as(x,"phylo4")
	}
	f.d <- paste("%.", digits, "g", sep = "") 
	format<-match.arg(format,choices=c("simmap1.1","newick"))
	x.reorder<-reorder(x, "postorder")
	description.vector<-rep(NA,nEdges(x.reorder))
	mrcaNode<-MRCA(x.reorder,taxa.vector)
	for (edgeIndex in 1:nEdges(x.reorder)) {
		currentEdge<-edges(x.reorder)[edgeIndex,2]
		currentNode<-getNode(x.reorder,currentEdge)
		state=0
		if (currentNode==mrcaNode) {
			state=-1	
		}
		else if (length(which(descendants(x.reorder,mrcaNode,type="all")==currentNode))==1) { #the current node is in the set of descendants of the mrca node
			state=1
		}
		simmapLabel=""
		if (rootNode(x.reorder)!=currentNode) { #don't do this for the root node
			if (state!=-1) {
				simmapLabel=paste(":{",state,",",sprintf(f.d, edgeLength(x.reorder,currentNode)),"}",sep="")
			}
			else {
				if ((change.position>0 && change.position<1) || !suppress.zero) { #so, if we don't have to worry about trimming parts with zero state
					simmapLabel=paste(":{1,",sprintf(f.d, (1-change.position)*edgeLength(x.reorder,currentNode)),":0,",sprintf(f.d, change.position*edgeLength(x.reorder,currentNode)),"}",sep="")
				}
				else {
					if(change.position==0) { #suppress.zero must be true to get this far
						simmapLabel=paste(":{1,",sprintf(f.d, edgeLength(x.reorder,currentNode)),"}",sep="")
					}	
					else { #change.position must be 1 and suppress.zero==true
						simmapLabel=paste(":{0,",sprintf(f.d, edgeLength(x.reorder,currentNode)),"}",sep="")
					}
				}
			}
			if (format=="newick") {
				simmapLabel=paste(":",sprintf(f.d, edgeLength(x.reorder,currentNode)),sep="") #in this case, override simmap. This is mostly for debugging
			}
		}
		if(nodeType(x.reorder)[getNode(x.reorder,currentEdge)]=="tip") {
			description.vector[edgeIndex]=paste(names(getNode(x.reorder,currentNode)),simmapLabel,sep="")
		}
		else { #internal node, perhaps even the root
			tmpDescription="("
			childrenNodes<-children(x.reorder,currentNode)
			for (childIndex in 1:length(childrenNodes)) {
				tmpDescription=paste(tmpDescription,description.vector[which(edges(x.reorder)[,2]==childrenNodes[childIndex])],sep="")
				if (childIndex<length(childrenNodes)) {
					tmpDescription=paste(tmpDescription,",",sep="")
				}
			}
			tmpDescription=paste(tmpDescription,")",simmapLabel,sep="")
			description.vector[edgeIndex]=tmpDescription
		}
	}
	return(description.vector[nEdges(x.reorder)]) #last element is the root
}