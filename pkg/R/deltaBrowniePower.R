streeBrlen<-function(n,type=c("balanced","left","right")) {
	kappa=1
	type<-match.arg(type)
	if(type=="balanced") {
		kappa=0	
	}
	return(kappaTree(compute.brlen(stree(n,type)),kappa))
}

blMultiplier<-function(phy,rate, change, shape){
	phy4<-as(phy, 'phylo4')
	#print(edgeLength(phy4))
	focal.node<-rootNode(phy4) #just to initialize it
	if(change=="root") {
		focal.node<-children(phy4, rootNode(phy4))[1]
	}
	else if (change=="cherry") {
		#find the focal node that has only 2 TERMINAL descendants
		
		#MAKE IT TAXON 1 ND SISTER
	}
	else if (change=="quarter") {
		#find the internal node that has ntax/4 TERMINAL descendants
		#MAKE IT SO THAT TAXON 1 IS PART OF THIS CLADE (ONLLY BALANCED)
	}
	focal.edge<-edgeLength(phy4)[which(edgeId(phy4)==getEdge(phy4,focal.node))]
	transformed.focal.edge<-0.5*(focal.edge + rate*focal.edge)
	edgeLength(phy4)[which(edgeId(phy4)==getEdge(phy4,focal.node))]<-transformed.focal.edge
	
	rest.clade<-descendants(phy4, focal.node, type="all")
	
	
	if(length(rest.clade)>1){
		for(node.index in 1:length(rest.clade)){
			new.focal.node=rest.clade[node.index]
			edgeLength(phy4)[which(edgeId(phy4)==getEdge(phy4,new.focal.node))]<-edgeLength(phy4)[which(edgeId(phy4)==getEdge(phy4,new.focal.node))]*rate
		}	
	}
	

	
	
#	print(edgeLength(phy4))
	return(as(phy4,"phylo"))	
}

ntax.vector=c(2^4, 2^5, 2^6, 2^7, 2^8)
shape.vector=c("balanced", "right")
rate.vector=c(0.1, 0.5, 1, 2, 10)
change.position.vector=c("root", "quarter", "cherry")
reps.per.combination=10


for(rep in 1:reps.per.combination){
	for(ntax.index in 1:length(ntax.vector)){
		ntax=ntax.vector[ntax.index]
	
		for(shape.index in 1:length(shape.vector)){
			shape=shape.vector[shape.index]
		
			for(rate.index in 1:length(rate.vector)){
				rate=rate.vector[rate.index]
				
				for(change.position.index in 1:length(change.position.vector)){
					change=change.position.vector[change.position.index]
					fileNameRoot<-paste("ntax",ntax,"shape",shape,"rate",rate,"change",change,"rep",rep,sep="",collapse="")
					batchFileName<-paste(fileNameRoot,".R",sep="",collapse="")
					cat("library(geiger)\n","library(RBrownie)\n",sep="",file=batchFileName,append=FALSE)
					cat("ntax=",ntax,"\nshape='",shape,"'\nrate=",rate,"\nchange='",change,"'\nrep=",rep,"\n",sep="",file=batchFileName,append=TRUE)
					cat("source('deltaBrownie.R')","\n",sep="",file=batchFileName,append=TRUE)
					cat("streeBrlen(ntax, type=shape)->phy\n",sep="",file=batchFileName,append=TRUE)
					
					
					#TO DO: KEEP ADDING COMMANDS FOR THE BATCH FILE 
					
					#MAKE SURE EACH DATA.FRAME IS SAVED AS AN R OBJECT IN ANOTHER FILE 
					
					#TO DO: SYSTEM(EZSUB R CMD BATCH BATCHFILENAME)
					
					
				}
			}
		}
	}	
}






