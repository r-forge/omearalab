#this will export the correct formats for cytoscape
source("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/UnifiedApproachScripts/V6_UtilityFns.R")
load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/all.results.cleaned.Rsave") 

plot.reduced=TRUE

#FOR PLOTTING REDUCED NETWORK, only take best values
if (plot.reduced) {
  highlevel.dataframe<-highlevel.dataframe[1:10,]
}



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
#qDataFrame$q=qDataFrame$q/min(qDataFrame$q)
print(qDataFrame)

if(plot.reduced) {
  #make it so only care about positions 1, 3, and 4. So only care about rates that differ here
  qDataFrame.pruned<-data.frame()
  for (q.row in sequence(dim(qDataFrame)[1])) {
     from.vector<-stringToVector(as.character(qDataFrame[q.row, ]$from))
     to.vector<-stringToVector(as.character(qDataFrame[q.row, ]$to))
     if(vectorMismatchExcludePositions(from.vector, to.vector, c(2,5,6))==1) {
       new.from=paste(from.vector[1],"x", from.vector[3], from.vector[4], "xx", sep="")
       new.to=paste(to.vector[1],"x", to.vector[3], to.vector[4], "xx", sep="")
       qDataFrame.pruned<-rbind(qDataFrame.pruned,data.frame(from=new.from, to=new.to, q=qDataFrame[q.row,]$q))

     }
  }
  qDataFrame.pruned<-unique(qDataFrame.pruned)
  qDataFrame<-qDataFrame.pruned
}

write.table(qDataFrame,file="/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/V7_Cytoscape_edges_NoBoring.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


#get info on number of combos
sourcetraits="/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv"
file<-sourcetraits
phy<-"/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/floral_1.nex"
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
  currentVector<-as.numeric(unlist(subdata[i,2:7]))
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
#print(ndf3)

if (plot.reduced) {
  ndf3$combo<-as.character(ndf3$combo)
  ndf4<-data.frame()
  for (node.row in sequence(dim(ndf3)[1])) {
    combo.vector<-strsplit(ndf3$combo[node.row],"")[[1]]
    new.combo<-paste(combo.vector[1], "x", combo.vector[3], combo.vector[4], "xx", sep="")
    if (node.row==1) {
      ndf4<-rbind(ndf4, ndf3[node.row, ])
      ndf4[1,]$combo<-new.combo
    }
    else {
      if(length(which(ndf4$combo == new.combo))==0) {
        ndf4<-rbind(ndf4, ndf3[node.row, ])
        ndf4$combo[dim(ndf4)[1]]<-new.combo
      }
      else {
        matching.row<-which(ndf4$combo == new.combo)
        ndf4$number_taxa_with_state[matching.row]<-ndf4$number_taxa_with_state[matching.row]+ndf3$number_taxa_with_state[node.row]
        ndf4$proportion_taxa_with_state[matching.row]<-ndf4$proportion_taxa_with_state[matching.row]+ndf3$proportion_taxa_with_state[node.row]
        
      }
    }
  }
  ndf4<-cbind(ndf4, time_until_transition=NA)
  for (node.row in sequence(dim(ndf4)[1])) {
    print(qDataFrame$q[which(qDataFrame$from == ndf4$combo[node.row])])
     ndf4$time_until_transition[node.row]<-1/sum(qDataFrame$q[which(qDataFrame$from == ndf4$combo[node.row])])
  }
  ndf4<-cbind(ndf4, extinction.fraction=ndf4$death/ndf4$birth)
  ndf3<-ndf4
}

write.table(ndf3,file="/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/V7_Cytoscape_nodes_NoBoring.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

#print(subdata)
