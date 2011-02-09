library(phylobase)
generate.simmap<-function(x, taxa.vector, change.position=0.5, digits=10, suppress.zero=FALSE) {
	f.d <- paste("%.", digits, "g", sep = "") 
	x.reorder<-reorder(x, "postorder")
	description.vector<-rep(NA,nEdges(x.reorder))
	mrcaNode<-MRCA(x.reorder,taxa.vector)
	nodeToEdgeIndex<-matrix(nrow=0,ncol=2)
	print(edges(x.reorder))
	for (edgeIndex in 1:nEdges(x.reorder)) {
		currentEdge<-edges(x.reorder)[edgeIndex,2]
		currentNode<-getNode(x.reorder,currentEdge)
		nodeToEdgeIndex<-rbind(nodeToEdgeIndex,c(currentNode,edgeIndex))
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
		}
		if(nodeType(x.reorder)[getNode(x.reorder,currentEdge)]=="tip") {
			description.vector[edgeIndex]=paste(names(getNode(x.reorder,currentNode)),simmapLabel,sep="")
		}
		else { #internal node, perhaps even the root
			tmpDescription="("
			print(nodeToEdgeIndex)
			childrenNodes<-children(x.reorder,currentNode)
			for (childIndex in 1:length(childrenNodes)) {
				#tmpDescription=paste(tmpDescription,description.vector[nodeToEdgeIndex[which(nodeToEdgeIndex[,1]==childrenNodes[childIndex])],2],sep="")
				if (childIndex<length(childrenNodes)) {
					tmpDescription=paste(tmpDescription,",",sep="")
				}
			}
			tmpDescription=paste(tmpDescription,")",simmapLabel,sep="")
			description.vector[edgeIndex]=tmpDescription
		}
	}
	print(description.vector)
}