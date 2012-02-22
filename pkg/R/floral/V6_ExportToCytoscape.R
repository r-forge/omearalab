#this will export the correct formats for cytoscape
source("~/Sites/RunsJan2012/UnifiedApproachScripts/V6_UtilityFns.R")
load("~/Sites/RunsJan2012/Summaries/Highlevel.dataframe.withrates.Rsave") #replace with frame of all data

bestValues<-highlevel.dataframe[which.max(highlevel.dataframe$AICweight),]
meanValues<-bestValues
for ( i in 11:length(meanValues) ) {
  meanValues[i]<-weighted.mean(highlevel.dataframe[,i],highlevel.dataframe$AICweight,na.rm=TRUE)
}

#first the edges
qValuesPositions<-grep("^q",names(meanValues))
qDataFrame<-data.frame()
for(i in 1:length(qValuesPositions)) {
  position<-qValuesPositions[i]
   if(length(grep("disallowed",names(meanValues)[position]))==0) {
    labels<-sub("q","",strsplit(names(meanValues)[position],"_")[[1]])
    rate<-meanValues[1,position]
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
write.table(qDataFrame,file="~/Sites/RunsJan2012/Summaries/Cytoscape_edges.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


#get info on number of combos
sourcetraits="../SourceData/Stebbins_prunenoper25i2012BCO.csv"
file<-sourcetraits
phy<-"../SourceData/floral_1.nex"
tree<-read.nexus(phy)
data<-read.csv(file)
colnamesVector<-colnames(data)
names<-colnamesVector[2:length(colnamesVector)]
#keep only 1 (1st 25%) and 2 (2nd 25%) in Selection column
#subdata<-subset(data,data$Selection<=2)
#already done, so just
subdata<-data #to minimize recoding

#delete taxa with missing data anywhere
for (colToExamine in 2:length(colnamesVector)) {
	subdata<-subdata[subdata[,colToExamine]!="?",]
}

comboCounts<-rep(0,2^S)
for(i in 1:dim(subdata)[1]) {
  currentVector<-as.numeric(unlist(subdata[i,2:8]))
  comboDecimal<-comboAsDecimal(currentVector,S)
  comboCounts[comboDecimal]<-comboCounts[comboDecimal]+1
}
comboProportions<-comboCounts/(sum(comboCounts))

#nodesDataFrame<-data.frame()
#now the nodes
try(for (comboDecimal in 1:2^S) {
  comboName<-comboAsBinaryString(comboDecimal,S)
  birthRate<-as.vector(meanValues[1,which(names(meanValues)==paste("lambda",comboName,sep=""))])
  deathRate<-as.vector(meanValues[1,which(names(meanValues)==paste("mu",comboName,sep=""))])
  diversificationRate<-birthRate-deathRate
  turnoverRate<-birthRate+deathRate
  tmp.dataFrame<-data.frame(comboName,birthRate,deathRate,diversificationRate,turnoverRate,comboCounts[comboDecimal],comboProportions[comboDecimal],vectorMismatch(comboAsBinaryVector(1,S),comboAsBinaryVector(comboDecimal,S)))
  names(tmp.dataFrame)<-c("combo","birth","death","diversification","turnover","number_taxa_with_state","proportion_taxa_with_state","steps_from_root")
 # tmp.dataFrame<-data.frame(5)
if (comboDecimal==1) {
   ndf2<-tmp.dataFrame
  }
  else {
    print(comboDecimal)
    print(tmp.dataFrame)
    ndf2<-rbind(ndf2,tmp.dataFrame)
    rm(tmp.dataFrame)
    print("done rbinding")
  }
  print(dim(ndf2))
  ndf3<-ndf2
  print("done if else")
})
print("done combo decimal")
print(ndf3)
write.table(ndf3,file="~/Sites/RunsJan2012/Summaries/Cytoscape_nodes.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

print(subdata)
