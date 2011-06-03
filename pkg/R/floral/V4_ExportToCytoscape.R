#this will export the correct formats for cytoscape
source("~/Sites/Floral/RunsApril2011/UnifiedApproachScripts/V4_UtilityFns.R")
load("~/Sites/Floral/RunsApril2011/Summaries/RateSummaryT 1 D 5 .Rsave") #replace with frame of all data
sdf<-summary.dataframe #just to save on typing
rm(summary.dataframe)
bestValues<-sdf[which.max(sdf$AICweight),]

#first the edges
qValuesPositions<-grep("^q",names(bestValues))
qDataFrame<-data.frame()
for(i in 1:length(qValuesPositions)) {
  position<-qValuesPositions[i]
   if(length(grep("disallowed",names(bestValues)[position]))==0) {
    labels<-sub("q","",strsplit(names(bestValues)[position],"_")[[1]])
    rate<-bestValues[1,position]
    names(rate)<-c("q")
    if(i==1) {
      qDataFrame<-data.frame(labels[1],labels[2],rate)
    }
    else {
      qDataFrame<-rbind(qDataFrame,data.frame(labels[1],labels[2],rate))
    }
   }
}
names(qDataFrame)<-c("from","to","q")
qDataFrame$q=qDataFrame$q/min(qDataFrame$q)
print(qDataFrame)
write.table(qDataFrame,file="~/Sites/Floral/RunsApril2011/Summaries/Cytoscape_edges.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

nodesDataFrame<-data.frame()
#now the nodes
for (comboDecimal in 1:S^7) {
  comboName<-comboAsBinaryString(comboDecimal,S)
  birthRate<-as.vector(bestValues[1,which(names(bestValues)==paste("lambda",comboName,sep=""))])
  deathRate<-as.vector(bestValues[1,which(names(bestValues)==paste("mu",comboName,sep=""))])
  diversificationRate<-birthRate-deathRate
  turnoverRate<-birthRate+deathRate
   tmp.dataFrame<-data.frame(comboName,birthRate,deathRate,diversificationRate,turnoverRate)
  names(tmp.dataFrame)<-c("combo","birth","death","diversification","turnover")
  tmp.dataFrame<-data.frame(5)
 if (comboDecimal==1) {
   nodesDataFrame<-tmp.dataFrame
  }
  else {
    nodesDataFrame<-rbind(nodesDataFrame,tmp.dataFrame)
  }
  print(nodesDataFrame)
}