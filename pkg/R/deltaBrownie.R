#this code should find the best break among alternate morphological models

#library(phylobase)
#library(RBrownie)
#library(geiger) 


#for three-taxa with real names example (short example)
#tip.label<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp")
#tree1<-rcoal(5, tip.label)
#tip.label2<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi")
#tree2<-rcoal(5, tip.label)


#tip.label<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp")
#completeTree<-rcoal(5, tip.label)

#completeBipesData<-read.table("/Users/halamillo/Desktop/completeBipes.txt", row.names=1)

#bipes.results<-iterateNonCensored(tree1, completeBipesData)

#tip.label2<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi")
#partialTree<-rcoal(4, tip.label2)

#partialBipesData<-read.table("/Users/halamillo/Desktop/partialBipes.txt", row.names=1)

#tip.label3<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp", "Bipes_sp2")
#extraTree<-rcoal(6, tip.label3)

#extraBipesData<-read.table("/Users/halamillo/Desktop/extraBipes.txt", row.names=1)


#data = runif((nTips(tree1)), 1, 10)
#data
#phy<-tree1
#bipes.results<-iterateNonCensored(phy, extraBipesData)

#for example with dipsadine tree and junk data
#phy<-read.tree("/Users/halamillo/Desktop/NonAcrochMonoRegCons2.phy")
#data<-read.table("/Users/halamillo/Desktop/dipsmorphdatatabdel.txt", row.names=1)

#this is a function to convert from phylo or phylo4 to a newick string, by default in simmap1.1 format. 
#By giving it a vector of taxon labels, it can give everything in the smallest clade containing those taxa
#a simmap state of 1, everything else a state of 0, and the edge leading to that clade can be broken up: 
#by default, half is in each state, but if change.position=0.75, the first 3/4 is in state 0 and the second 1/4
#is in state 1.
#If that edge is set to be entirely in one state (change.position==1 or change.position==0), there will be 
#one state with zero length on that branch. This can be optionally deleted using suppress.zero=TRUE

generate.simmap<-function(x, taxa.vector, change.position=0.5, digits=10, suppress.zero=FALSE,format="simmap1.1") {
	if(class(x)!="phylo4") {
		x<-as(x,"phylo4")
	}
	f.d <- paste("%.", digits, "g", sep = "") 
	format<-match.arg(format,choices=c("simmap1.1","newick"))
	x.reorder<-reorder(x, "postorder")
	description.vector<-rep(NA,nEdges(x.reorder))
	mrcaNode<-MRCA(x.reorder,taxa.vector)
	node_number<-mrcaNode
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


tree<-description.vector[nEdges(x.reorder)]
names(tree)<-c(paste(node_number, sep=""))
return(tree) #last element is the root

}



iterateNonCensored<-function (phy, data, name.check=TRUE) {
	
	if(name.check){
		phy<-as(phy,"phylo")
		
#check that taxa match between data and tree
		checked.object<-name.check(phy,data)
		
		if (checked.object[1] == "OK") { #pillaged kindly from Banbury
			cat("Phylogeny and character matrix are in agreement...moving-on with the rate analysis!")
			
			
			
		} else {
			warning("Phylogeny and character matrix are mismatched:\n")
			cat("missmatched taxa have been dropped from analysis\n")
#cat("(In case you are curious check your working directoy for files(s) with data that didn't match)\n")
			
			if (length(checked.object$Tree.not.data) > 0) {  
#a<-checked.object$Tree.not.data
#write(a, file="Tree.not.data.txt")
				phy<-drop.tip(phy, checked.object$Tree.not.data)
			}
			
			if (length(checked.object$Data.not.tree) > 0) {
#b<-checked.object$Data.not.tree
#write(b, file="Data.not.tree.txt")
				which(rownames(data) %in% checked.object$Data.not.tree)->rows
				as.data.frame(data[-rows,])->data1
				colnames(data1)<-colnames(data)
				rownames(data1)<-rownames(data)[-rows]
				data<-data1
			}
		}
	}
	phy<-as(phy, "phylo4d") #turn ape phylogeny into phylo4d class
	browniePhy<-as(phy, "brownie") #turn phy into a brownie object
	brownie.tree1<-addData(browniePhy, tip.data=data, dataTypes=contData())  #add data to browniePhy
	brownie.tree1@order<-"unknown"  #patch so descendants() for the root node works
###for when we want to index all clade combos###
	Nodes1<-grep("tip",nodeType(brownie.tree1))  #pulls out the nodes that are internal
	Nodes2<-grep("internal",nodeType(brownie.tree1))  #pulls out the nodes that are internal
	Nodes<-c(Nodes1, Nodes2)
	allDec<-vector("list", length(Nodes)+1)  #creates list with as many places as there are Nodes plus 1 to make space for the TAXSET_all taxon set
	names(allDec)<-paste("TAXSET_", Nodes,sep="") #all names are added the "TAXSET_" necessary for the brownie object to recognize them
	
	for (ii in 1:length(allDec)){		#this for-loop creates the required (by RBrownie) TAXSET_all taxset; thank you Conrad for patching the rootNode mess!
		if (is.na(names(allDec[ii]))) {
			names(allDec)[ii]<-paste("TAXSET_all")
			allDec[[ii]]<-descendants(brownie.tree1, rootNode(brownie.tree1), type=c("tips"))
			taxasets(brownie.tree1, taxnames=names(allDec[ii]))<-labels(allDec[[ii]])
		}
	}
	
	for (i in 1:length(Nodes)) {     #this for-loop creates the required (by RBrownie) TAXSETS_ for each inner node of the tree
		allDec[[i]]<-descendants(brownie.tree1, Nodes[i], type=c("tips"))
		taxasets(brownie.tree1, taxnames=names(allDec[i]))<-labels(allDec[[i]])
	}
	
	nombres<-vector("list", length(allDec))
	
	for (i in 1:length(allDec)){  #pulls out taxa that will be used in the manufacutre of taxa.vectors
		nombres[[i]]<-labels(allDec[i][[1]])
	}
	for (i in 1:length(nombres)){ # create the empty list of taxa.vectors
		names(nombres)[i]<-paste("taxa.vector.",i, sep="")
	}
	cat("nombres")
####so now here use each of those taxa.vectors to generate the other simmap phylogenies
#### and read in each simmap formatted tree into an object
	
	trees.list<-c()
	trees.list.phy<-vector("list", length(nombres))
	for (i in 1:length(nombres)){
		trees.list[i]<-generate.simmap(phy, nombres[i][[1]])
		trees.list.phy[i]<-read.simmap(text=trees.list[i])
	}
	cat("trees.list  ")	
####turn each tree into a phylo4d_ext class
	phy.ext.list<-vector("list", length(trees.list.phy))
	for(i in 1:length(trees.list.phy)){
		phy.ext.list[i]<-phyext(trees.list.phy[[i]])		}	
	cat("phy.ext.list  ")
####turn each tree into a brownie class	
	phy.brownie.list<-vector("list", length(phy.ext.list))
	for(i in 1:length(phy.ext.list)){
		phy.brownie.list[i]<-brownie(phy.ext.list[i])		}
	cat("phy.brownie.list  ")	
#####add the data to each tree
	phy.brownie.list.w.data<-vector("list", length(phy.brownie.list))
	for(i in 1:length(phy.brownie.list)){
		phy.brownie.list.w.data[i]<-addData(phy.brownie.list[i], tip.data=data, dataTypes=contData())
	}
	cat("phy.brownie.list.w.data  ")
	
#generate the names list of where the tree brakes
	NodeShift=c()
	
	for (i in 1:length(nombres)){
		NodeShift[i]<-names(generate.simmap(phy, nombres[i][[1]]))
		NodeShift[i+1]<-c("NA")
	}
	
####now generate the "all" taxaset
	all_taxa<-grep("[a-z]", tipLabels(phy.brownie.list.w.data[[1]]), value=TRUE)  ##NOTE this is only for the junk tree -- need to fix this for the real species names 
	for(i in 1:length(phy.brownie.list.w.data)){
		taxasets(phy.brownie.list.w.data[[i]], taxnames="all")<-all_taxa		}
	cat("all_taxa.grep  ")
	all.test.results1<-data.frame() #create an empty data frame as a repository of the final results
	all.test.results2<-data.frame() #create an empty data frame as a repository of the final results
	
	
	
	all.test.results1<-runNonCensored(phy.brownie.list.w.data, models=brownie.models.continuous()[2], treeloop=T, charloop=T)
	all.test.results2<-runNonCensored(phy.brownie.list.w.data[[1]], models=brownie.models.continuous()[1], treeloop=T, taxset="all")
	
	
	
	newRow<-data.frame(Tree=all.test.results2$Tree, Tree.weight =all.test.results2$"Tree weight", Tree.name=all.test.results2$"Tree name", Char=all.test.results2$Char, Model=all.test.results2$Model, LnL=all.test.results2$"-LnL", AIC=all.test.results2$AIC, AICc=all.test.results2$AICc, AncState=all.test.results2$AncState, Rate_in_state_0=all.test.results2$BMrate, Rate_in_state_1=all.test.results2$BMrate)
	
	colnames(newRow)=c("Tree", "Tree weight", "Tree name", "Char", "Model", "-LnL", "AIC", "AICc", "AncState", "Rate_in_state_0", "Rate_in_state_1")
	
	all.test.results<-rbind(all.test.results1, newRow)
	
	all.test.results<-cbind(all.test.results,NodeShift)
	
	
	dAIC=all.test.results$AIC-min(all.test.results$AIC)
	dAICc=all.test.results$AICc-min(all.test.results$AICc)
	
	AICweightraw=exp(-0.5*dAIC)
	AICweight=AICweightraw/sum(AICweightraw)
	
	
	AICcweightraw=exp(-0.5*dAICc)
	AICcweight=AICcweightraw/sum(AICcweightraw)
	
	all.test.results<-cbind(all.test.results, dAIC, AICweight)
	all.test.results<-cbind(all.test.results,dAICc, AICcweight)
	
	colnames(all.test.results)[6]<-"-LnL"
	
	result.object<-vector("list",3)
	names(result.object)<-c("all.test.results", "rate.at.shift", "tree.with.shift")
	
	result.object$all.test.results<-all.test.results
	
	result.object$rate.at.shift<-result.object$all.test.results[which(result.object$all.test.results$AICcweight==max(result.object$all.test.results$AICcweight)),]
	
	result.object$tree.with.shift<-phy.brownie.list.w.data[[result.object$all.test.results$Tree[which((result.object$all.test.results$AICcweight)==max(result.object$all.test.results$AICcweight))]]]
#result.object<-c(all.test.results, trees.list.phy)
#return(all.test.results)
	return(result.object)
	
}




#######Function to generate a balanced, left-leaning, or right-leaning phylogeny with branch lengths
streeBrlen<-function(n,type=c("balanced","left","right")) {
	kappa=1
	type<-match.arg(type)
	if(type=="balanced") {
		kappa=0	
	}
	return(kappaTree(compute.brlen(stree(n,type)),kappa))
}



######Function to multiply the branch lenghts of a clade that can happen at the root, tipward (i.e. cherry), or that comprises 1/4 of the species
blMultiplier<-function(phy,rate, change, shape){
	phy4<-as(phy, 'phylo4')
#print(edgeLength(phy4))
	focal.node<-rootNode(phy4) #just to initialize it
	
	if(change=="root") {
		focal.node<-children(phy4, rootNode(phy4))[1]
	}
	else if (change=="cherry") {
#find the focal node that has only 2 TERMINAL descendants
		
		if(shape=="right"||shape=="left"){
			for(i in 1:nNodes(phy4)){
				while(length(descendants(phy4, focal.node, type=c("tips")))>2){
					focal.node<-descendants(phy4, focal.node, type=c("children"))[which(descendants(phy4, focal.node, type=c("children")) %in% as.numeric(labels(grep("internal", nodeType(phy4)[children(phy4, focal.node)], value=TRUE))))]
				}
			}
		}
		else if(shape=="balanced"){
			focal.node<-ancestor(phy4, descendants(phy4, focal.node, type=c("tips"))[1])
		}
		
	}
	else if (change=="quarter") {
#find the internal node that has ntax/2 TERMINAL descendants
#MAKE IT SO THAT TAXON 1 IS PART OF THIS CLADE (ONLLY BALANCED)
		quarter.specs<-((nTips(phy4))/4)
		desc.down.node<-c()	
		
		if(shape=="right"||shape=="left"){
			for(i in 1:nNodes(phy4)){
				while(length(descendants(phy4, focal.node, type=c("tips")))>2){
					focal.node<-descendants(phy4, focal.node, type=c("children"))[which(descendants(phy4, focal.node, type=c("children")) %in% as.numeric(labels(grep("internal", nodeType(phy4)[children(phy4, focal.node)], value=TRUE))))]
				}
			}
		}
		
		else if(shape=="balanced"){	
			while(length(descendants(phy4, focal.node, type=c("tips")))>1){
				focal.node<-descendants(phy4, focal.node, type=c("children"))[1]
			}
		}
		
		while(length(desc.down.node)<quarter.specs){
			sibs.focal.node<-siblings(phy4, focal.node)
			f.n.and.sibs<-c(focal.node, sibs.focal.node)
			focal.node<-MRCA(phy4, f.n.and.sibs)
			desc.down.node<-descendants(phy4, focal.node, type=c("tips"))
		}
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