source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V6_UtilityFns.R")
library(diversitree)
library(sfsmisc) 
library(ape)
library(partitions) #for converting from binary back to decimal
library(corpcor)
load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/all.results.cleaned.Rsave") 

bestValues<-highlevel.dataframe[which.max(highlevel.dataframe$AICweight),]
meanValues<-bestValues
for ( i in 11:length(meanValues) ) {
  meanValues[i]<-weighted.mean(highlevel.dataframe[,i],highlevel.dataframe$AICweight,na.rm=TRUE)
}


sourcetraits="/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv"
file<-sourcetraits
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
names(comboCounts)<-sapply(c(1:64),comboAsBinaryString, S=S)

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

colfunc <- colorRampPalette(c("blue", "red"))(100)


plot(c(-5,64), c(-.2*max(comboProportions), max(comboProportions)), bty="n", xlab="", ylab="", xaxt="n", yaxt="n", type="n")
lines(x=c(-5,64), c(0,0), col="gray")
state.names<-c("Perianth", "Fusion", "Symmetry", "Stamens", "Syncarpy", "Ovary")
text(x=rep(0,6),y=seq(from=-1, to=-6, length.out=6)*.2*max(comboProportions)/6, labels=state.names, pos=2, cex=0.5)
state.id<-c(1:64)
state.id<-state.id[order(comboProportions,decreasing=TRUE)]
for (i in 1:64) {
  state<-state.id[i]
  col.value<-round(((ndf3$diversification[state]-min(ndf3$diversification))/(max(ndf3$diversification)-min(ndf3$diversification))),2)
  col<-colfunc[100*col.value]
  if(col.value==0) {
    col="blue" 
  }
  print(c(col.value, col))
 lines(x=rep(i,2), y=c(0, comboProportions[state]),col=col) 
 states<-comboAsBinaryVector(state, S)
 for (j in sequence(length(states))) {
    if(states[j]==0) {
      points(i,j*(-.2*max(comboProportions))/6, pch=".", bg=col, col=col) 
    }
 }
}

#axis(side=2, pos=0, at=seq(from=0,to=max(comboProportions), length.out=3), labels=c("","",round(max(comboProportions),2)))

save(list=ls(), file="ComponentsForPlottingHistogram.RData", compress=TRUE)