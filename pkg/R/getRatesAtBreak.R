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
		else { 
      rate.choice[each.node,1:2]=1 #ratechoice.pre where row is each.node, where column 1 and column 2 = rate.choice
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


#result.iterateNonCensored<-bipes.results


#Stuff that is R but iscommented out is to fix
getModelAveragedRates<-function(phy, result.iterateNonCensored, change.position=0.5) {
	if(class(phy)!="phylo4") {
		phy<-as(phy,"phylo4")
	}
	node.vector<-c()
	nodeList<-nodeId(phy,type=c("all"))	
	rate.matrix<-matrix(data=0, nrow=length(nodeList), ncol=2)  #this is where we will stick the model-averaged rates
	for(modelIndex in 1:length(result.iterateNonCensored$all.test.results$NodeShift)) { #this goes by row and NodeShift does not match the order of the rows see model.rates[3]	
		model.rates<-rep(NA,2)
		model.rates[1]<-result.iterateNonCensored$all.test.results$Rate_in_state_0[modelIndex] #CHECK THIS SYNTAX, make sure it works for BM1, too
		model.rates[2]<-result.iterateNonCensored$all.test.results$Rate_in_state_1[modelIndex] #ditto
		AICcweight<-result.iterateNonCensored$all.test.results$AICcweight[modelIndex]
		BMS.node<-result.iterateNonCensored$all.test.results$NodeShift[modelIndex] 
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


avgRateNodeShiftInfo<-function(avg.rates.object){
		rates.table<-avg.rates.object		
		rate.diff.matrix<-matrix(data=0, nrow=dim(rates.table)[1], ncol=1)		
		for(rate.index in 1:dim(rates.table)[1]){
			rate.diff.matrix[rate.index]<-(rates.table[rate.index, 1] - rates.table[rate.index, 2])
		}	
		rate.diff.matrix<-abs(rate.diff.matrix)
		Avg.rt.shift<-which(rate.diff.matrix == max(rate.diff.matrix))
		best.rates<-rates.table[Avg.rt.shift, 1:2]
		proportional.change<-max(rates.table[Avg.rt.shift, 1:2])/min(rates.table[Avg.rt.shift, 1:2])
		which(max(rates.table[Avg.rt.shift, 1:2]) == rates.table[Avg.rt.shift, 1:2]) ->pos
		if (pos == 2) {	
      type<-c("increase")
		}
		else{ 
      type<-("decrease")
		}
		result.object<-vector("list",4)
		names(result.object)<-c("avg.rate.node.shift", "best.avg.rates", "proportional.change", "type.of.change")
		result.object$avg.rate.node.shift<-Avg.rt.shift
		result.object$best.avg.rates<-best.rates
		result.object$proportional.change<-proportional.change
		result.object$type.of.change<-type
		return(result.object)
}

nodesToEdges<-function(phy) {
  if(class(phy)!="phylo4") {
    phy<-as(phy,"phylo4") 
  }
  nodeVector<-nodeId(phy)
  edgeVector<-rep(NA,length(nodeVector))
  edgesMatrix<-edges(phy)
  for (nodeIndex in 1:length(nodeVector)) {
     node<-nodeVector[nodeIndex]
     edgeVector[nodeIndex]<-which(edgesMatrix[,2]==node)
  }
  finaldataframe<-data.frame(nodes=nodeVector,edges=edgeVector)
  return(finaldataframe)
}

#this function just plots the lines and labels. It is normally called by plotAvgRates for one tree, but you can call
#  it separately with a nonzero yOffset to plot multiple trees in one plot

#todo: makes sure change.position is changed to match what it is in setting the break point in the multiple rate testing loop.Basically, if the split is set to happen a third of the way up the branch, allow for this (change getRatesAtBreak, too, to avoid assuming this happens at 0.5)
linesAvgRates<-function(phy,avgRates,globalMin=min(avgRates),globalMax=max(avgRates),yOffset=0,show.tip.label=TRUE,change.position=0.5,grayMax=0.8,lwd=2,textThreshold=1.5,rounddigits=1) {
  if(class(phy)!="phylo4") {
    phy<-as(phy,"phylo4") 
  }
  rateRatios01<-avgRates[,1]/avgRates[,2]
  rateRatios10<-avgRates[,2]/avgRates[,1]
  
  
  xxyy<-phyloXXYY(phy)
  nodesToEdges.df<-nodesToEdges(phy)
  vert0<-xxyy$segs$v0x
  
  x0 = xxyy$segs$v0x #modified from treePlot.R in phylobase
  y0 = xxyy$segs$v0y+yOffset
  x1 = xxyy$segs$h1x
  y1 = xxyy$segs$h1y+yOffset
 # y1=y0
  xmid<-change.position*x0+(1-change.position)*x1
  ymid<-change.position*y0+(1-change.position)*y1
  for (nodeIndex in 1:dim(avgRates)[1]) {
    node<-as.numeric(row.names(avgRates)[nodeIndex])
    #print(xxyy$eorder)
    #print(nodesToEdges.df[,1])
    #print(node)
    segmentPointer<-which(xxyy$eorder==nodesToEdges.df[(which(nodesToEdges.df[,1]==node)),2])
    rate0<-avgRates[nodeIndex,1]
    rate0proportion<-(rate0-globalMin)/(globalMax-globalMin)
    rate1<-avgRates[nodeIndex,2]
    rate1proportion<-(rate1-globalMin)/(globalMax-globalMin)
    rate0rescaled<-rate0/min(rate0,rate1)
    rate1rescaled<-rate1/min(rate0,rate1)
    #remember, with gray, 0=black, 1=white. 
    #vertical lines
    lines(x=c(x0[segmentPointer],x0[segmentPointer]),y=c(y0[segmentPointer],y1[segmentPointer]),col=gray((1-rate0proportion)*grayMax),lwd=lwd)
    #horizontal lines
    lines(x=c(x0[segmentPointer],xmid[segmentPointer]),y=c(y1[segmentPointer],y1[segmentPointer]),col=gray((1-rate0proportion)*grayMax),lwd=lwd)
    lines(x=c(xmid[segmentPointer],x1[segmentPointer]),y=c(y1[segmentPointer],y1[segmentPointer]),col=gray((1-rate1proportion)*grayMax),lwd=lwd)
    writeText<-FALSE
    if (max(rate1rescaled,rate0rescaled)>textThreshold) {
      writeText<-TRUE 
    }
    if (writeText) {
      text(x=xmid[segmentPointer],y=y1[segmentPointer],labels=paste(round(rate0rescaled,digits=rounddigits),":",round(rate1rescaled,digits=rounddigits),sep=""),pos=3)
      rect(2,2,3,4, col=rgb(0,0,0,0.1), border=FALSE) #xleft, ybottom, xright, ytop; co=rgb(r,g,b,transparency)
       
    }
  }
}

plotAvgRates<-function(phy,avgRates,globalMin=min(avgRates),globalMax=max(avgRates),yOffset=0,show.tip.label=TRUE) {
  xxyy<-phyloXXYY(phy)
  plot(x=xxyy$xx,y=xxyy$yy,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  linesAvgRates(phy,avgRates,globalMin,globalMax,yOffset,show.tip.label)
}