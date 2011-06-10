#tip.label<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp")
#tree1<-rcoal(5, tip.label)> 
#data = runif((nTips(tree1)), 1, 10)
#data<-as.data.frame(data)
#rownames(data)<-tip.label
#bipes.results<-iterateNonCensored(tree1, data)
#bipes.results$tree.with.shift->phyext_tree
#plot(phyext_tree)

getRatesAtBreak<-function(phy, BMS.node) {

	node.vector<-c()
	#rate.choice<-c()
	nodeList<-nodeId(phy,type=c("all"))
	
	rate.choice<-matrix(data=NA, nrow=length(nodeList), ncol=2) 

	for(nodeIndex in 1:length(nodeList)){ 
		each.node<-nodeList[nodeIndex]
		if (length(which(ancestors(phy, each.node)==BMS.node))==0) {
			rate.choice[each.node,1:2]=0  #ratechoice.pre where row is each.node, where column 1 and column 2 = rate.choice
		}		
			else { rate.choice[each.node,1:2]=1 #ratechoice.pre where row is each.node, where column 1 and column 2 = rate.choice
			}

		if (each.node == BMS.node) {
	
			rate.choice[each.node,1:2]=c(0,1) #ratechoice.pre where row is each.node, where column 1 and column 2 = rate.choice    #but here it should be a 0 or a 1; decide which 
	
		}
	colnames(rate.choice)<-c("pre", "post")	
	cat(rate.choice)
	#rate.choice.pre<-append(rate.choice.pre, rate.choice, after=length(rate.choice.pre))
		
	}	
return(rate.choice)	
}	
