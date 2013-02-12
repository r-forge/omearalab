source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V6_UtilityFns.R")
library(diversitree)
library(sfsmisc) 
library(ape)
library(partitions) #for converting from binary back to decimal
library(igraph)
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

g <- graph.empty()
#nodesDataFrame<-data.frame()
#now the nodes
colors.rate<-c()
colors.anc<-c()
colors.focal<-c()
colors.edges<-c()
try(for (comboDecimal in 1:2^S) {
  comboName<-comboAsBinaryString(comboDecimal,S)
  g<-g+vertices(comboName)
  
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
   # print(comboDecimal)
  #  print(tmp.dataFrame)
    ndf2<-rbind(ndf2,tmp.dataFrame)
    rm(tmp.dataFrame)
#    print("done rbinding")
  }
#  print(dim(ndf2))
  ndf3<-ndf2
#  print("done if else")
})
#print("done combo decimal")

colfunc <- colorRampPalette(c("blue", "red"))(100)


#plot(c(-5,64), c(-.2*max(comboProportions), max(comboProportions)), bty="n", xlab="", ylab="", xaxt="n", yaxt="n", type="n")
#lines(x=c(-5,64), c(0,0), col="gray")
state.names<-c("Perianth", "Fusion", "Symmetry", "Stamens", "Syncarpy", "Ovary")
#text(x=rep(0,6),y=seq(from=-1, to=-6, length.out=6)*.2*max(comboProportions)/6, labels=state.names, pos=2, cex=0.5)
state.id<-c(1:64)
state.id<-state.id[order(comboProportions,decreasing=TRUE)]

#optimal focal is 0x11xx, so match that

for (i in 1:64) {
  state<-state.id[i]
  col.value<-round(((ndf3$diversification[state]-min(ndf3$diversification))/(max(ndf3$diversification)-min(ndf3$diversification))),2)
  col<-colfunc[100*col.value]
  if(col.value==0) {
    col="blue" 
  }
  colors.rate<-append(colors.rate, col)
  if(i==1) {
    colors.anc<-append(colors.anc, "red") 
  } else {
    colors.anc<-append(colors.anc, "black")
  }
  if (comboAsBinaryVector(i, S)[1]==0 && comboAsBinaryVector(i, S)[3]==1 && comboAsBinaryVector(i, S)[4]==1) {
    colors.focal<-append(colors.focal, "orange")
  } else {
    colors.focal<-append(colors.focal, "purple")
  }
  #print(c(col.value, col))
  #lines(x=rep(i,2), y=c(0, comboProportions[state]),col=col) 
  states<-comboAsBinaryVector(state, S)
  for (j in sequence(length(states))) {
    if(states[j]==0) {
      #points(i,j*(-.2*max(comboProportions))/6, pch=".", bg=col, col=col) 
    }
  }
}

weights<-c()
colors.edges<-c()
for (i in 1:2^S) {
  for (j in 1:2^S) {
    if(vectorMismatch(comboAsBinaryVector(i,S), comboAsBinaryVector(j,S))==1 & i>j) {
      g<-g + edge(i, j)
      weight<-3
      col<-"darkorchid1"
      if (colors.focal[i]=="orange" & colors.focal[j]=="orange") {
        col="orange" 
      }
      if (colors.focal[i] != colors.focal[j]) {
        col="aquamarine4"
        weight<-.1
      }
      colors.edges<-append(colors.edges, col)
      weights<-append(weights, weight)      
    }
  }
}

start<-matrix(nrow=2^S, ncol=2) 
for (i in 1:2^S) {
  if (colors.focal[i]=="orange") {
    start[i,]<-rnorm(2, -5, 0.001) 
  } else {
    start[i, ]<-rnorm(2, 5, 0.001) 
  }
}
colors.edges.gray<-colors.edges
colors.edges.gray[which(colors.edges.gray=="darkorchid1")]<-"lightgray"

colors.edges.grayscale<-colors.edges
colors.edges.grayscale[which(colors.edges=="darkorchid1")]<-"lightgray"
colors.edges.grayscale[which(colors.edges=="aquamarine4")]<-"darkgray"
colors.edges.grayscale[which(colors.edges=="orange")]<-"black"


size.focal<-rep(7, 2^S)
size.focal[which(colors.focal=="orange")]<-15

colors.focal.bw<-colors.focal
colors.focal.bw[which(colors.focal=="orange")]<-"white"
colors.focal.bw[which(colors.focal!="orange")]<-"black"


par(mfcol=c(1,2))
plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges, vertex.size=7)
plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges, vertex.size=15*comboProportions/max(comboProportions))

#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal.bw, vertex.label=NA, edge.arrow.size=0, edge.color="darkgray", vertex.size=7)


#plot(g, layout=layout.kamada.kawai(g, start=start),  vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges, vertex.size=7)
#plot(g, layout=layout.kamada.kawai(g, start=start),  vertex.color=colors.focal.bw, vertex.label=NA, edge.arrow.size=0, edge.color="darkgray", vertex.size=7)

#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges.grayscale, vertex.size=15*comboProportions/max(comboProportions))
#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal.bw, vertex.label=NA, edge.arrow.size=0, edge.color="darkgray", vertex.size=15*comboProportions/max(comboProportions))

#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.rate, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges.grayscale, vertex.size=7)
#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.rate, vertex.label=NA, edge.arrow.size=0, edge.color="darkgray", vertex.size=15*comboProportions/max(comboProportions))


#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights/weights), vertex.color="black", vertex.label=NA, edge.arrow.size=0, edge.color="darkgray")
#plot(g, layout=layout.circle(g),  vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges)
#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges.gray, vertex.size=7)
#plot(g, layout=layout.fruchterman.reingold(g, start=start, weights=weights), vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges, vertex.size=15*comboProportions/max(comboProportions))
#plot(g, layout=layout.kamada.kawai(g, start=start),  vertex.color=colors.focal, vertex.label=NA, edge.arrow.size=0, edge.color=colors.edges, vertex.size=15*comboProportions/max(comboProportions))
#plot(g, layout=layout.kamada.kawai(g, start=start),  vertex.color=colors.rate, vertex.label=NA, edge.arrow.size=0, edge.color="darkgray", vertex.size=7)
#plot(g, layout=layout.kamada.kawai(g, start=start),  vertex.color=colors.rate, vertex.label=NA, edge.arrow.size=0, edge.color="darkgray", vertex.size=size.focal)


