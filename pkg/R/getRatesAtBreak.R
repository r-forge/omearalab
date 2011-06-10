#tip.label<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp")
#tree1<-rcoal(5, tip.label)
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
	#cat(rate.choice)
	#rate.choice.pre<-append(rate.choice.pre, rate.choice, after=length(rate.choice.pre))
		
	}	
	return(rate.choice)	
}	


result.iterateNonCensored<-bipes.results


#Stuff that is R but iscommented out is to fix
getModelAveragedRates<-function(phy, result.iterateNonCensored) {
	node.vector<-c()
	nodeList<-nodeId(phy,type=c("all"))	
	rate.matrix<-matrix(data=0, nrow=length(nodeList), ncol=2)  #this is where we will stick the model-averaged rates

	for(modelIndex in 1:length(result.iterateNonCensored$all.test.results$NodeShift)) { #this goes by row and NodeShift does not match the order of the rows see model.rates[3]
		
		model.rates<-rep(NA,2)
		
		model.rates[1]<-result.iterateNonCensored$all.test.results$Rate_in_state_0[modelIndex] #CHECK THIS SYNTAX, make sure it works for BM1, too
		
		model.rates[2]<-result.iterateNonCensored$all.test.results$Rate_in_state_1[modelIndex] #ditto
		
		#model.rates[3]<-result.iterateNonCensored$all.test.results$NodeShift[modelIndex] #gives back value of node shift not the row 


		AICcweight<-result.iterateNonCensored$all.test.results$AICcweight[modelIndex]
		
		
		BMS.node<-result.iterateNonCensored$all.test.results$NodeShift[modelIndex] #does this need to be the NodeShift column?
		#if(result.iterateNonCensored$all.test.results$NodeShift[modelIndex]!="NA") {
		#	BMS.node<-result.iterateNonCensored$all.test.results$NodeShift[modelIndex]
		#}	
		
		rate.choice<-matrix(data=NA, nrow=length(nodeList), ncol=2)
		if(result.iterateNonCensored$all.test.results$Model[modelIndex]=="BM1") {
			rate.choice<-matrix(data=0, nrow=length(nodeList), ncol=2)
		}
		else {
			rate.choice<-getRatesAtBreak(phy,BMS.node)
		}
		
		for(nodeIndex in 1:length(nodeList)){
			each.node<-nodeList[nodeIndex]
			for (position in 1:2) {
				rate.matrix[nodeIndex,position] <- rate.matrix[nodeIndex,position] + ( AICcweight * model.rates[1+ rate.choice[nodeIndex,position] ] )
			}
		}
	}
	rownames(rate.matrix)<-nodeList
	return(rate.matrix)
}
