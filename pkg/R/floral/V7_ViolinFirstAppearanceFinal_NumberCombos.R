library(RColorBrewer)
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_NewSimulator.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_UtilityFns.R")
#dirs<-system("ls -1 | grep -v RData | grep -v '\\.R' | grep -iv pdf",intern=TRUE)
dirs<-c("just.the.50.final.runs")
#mypalette<-brewer.pal(8,"Dark2")
mypalette<-heat.colors(12,alpha=0.3)
original.data<-read.csv("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv",stringsAsFactors=FALSE)
observed.ntax<-dim(original.data)[1]
all.histories.list<-list()
combo.names<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")

GetFirstAppearance<-function(counts) {
  age<-NA
  nonzeros<-which(counts!=0)
  if(length(nonzeros)>0) {
    age<-nonzeros[1]-1 
  }
  if(age>136) {
  	age<-NA
  }
  return(age)
}

CountNonzeros<-function(x) {
	return(sum(x>0))
}


GetAgesOfOneThroughEightCombos<-function(history) {
	ages<-rep(NA,8)
	sums<-apply(history[,2:9],1, CountNonzeros)
	for (i in sequence(8)) {
		matches<-which(sums==i)
		if(length(matches)>0) {
			ages[i]<-matches[1]-1
			if(ages[i]>136) {
				ages[i]<-NA	
			}	
		}	
		
	}
	print(ages)
	return(ages)
}


observed.proportions<-rep(NA, length(combo.names))
for (i in sequence(length(combo.names))) {
  chars<-strsplit(combo.names[i],"")[[1]]
  small.dataset<-original.data[which(original.data[,1+1]==chars[1]),] #look at the relevant chars
  small.dataset<-small.dataset[which(small.dataset[,3+1]==chars[3]),]
  small.dataset<-small.dataset[which(small.dataset[,4+1]==chars[4]),]
  observed.proportions[i]<-dim(small.dataset)[1]/dim(original.data)[1]
}
names(observed.proportions)<-combo.names

summary.df<-cbind(data.frame(model="empirical", nrun=1, ntax.025=250000, ntax.med=250000, ntax.975=250000, dist.025=0, dist.med=0, dist.975=0), data.frame(t(round(observed.proportions, 2))))


subsample.proportion<-function(p,observed.ntax) {
  recovered.p<-rbinom(1,observed.ntax,p)/observed.ntax
  return(recovered.p)
}
par(mfcol=c(1,1))
plot(x=c(-10,150), y=c(0,9), yaxt="n", xlab="MY from root", bty="n", type="n",ylab="",xaxt="n")
axis(side=1, at=c(0, 136/2, 136))
for (dir.index in sequence(length(dirs))) {
  #print(dirs[dir.index])
  setwd(dirs[dir.index])
  ntax.vector<-c()
    file.list<-system("ls -1 | grep total151.RSave",intern=TRUE)
#  file.list<-system("ls -1 | grep total137.RSave",intern=TRUE)
  if(grepl("old", dirs[dir.index])) {
	file.list<-system("ls -1 | grep RSave", intern=TRUE)
}
  if (length(file.list)>0) {
    file.list<-paste("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/", dirs[dir.index], "/", file.list, sep="")
  }
  #  if(!grepl("stacey", dirs[dir.index])) {
  #    setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINALAdditional")
  #    setwd(dirs[dir.index])
  #    file.list2<-system("ls -1 | grep total137.RSave",intern=TRUE)
  #    if (length(file.list2)>0) {
  #      file.list2<-paste("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINALAdditional/", dirs[dir.index], "/", file.list2, sep="")
  #    }
  #    file.list<-c(file.list, file.list2)
  #  }
  # print(file.list)
  setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL")
  setwd(dirs[dir.index])
  distances<-c()
  for(combo.index in rev(sequence(8))) {
  		ages.of.first<-c()
  		for (i in sequence(length(file.list))) {
    	load(file.list[i])
    	ages.of.first<-append(ages.of.first, GetAgesOfOneThroughEightCombos(history)[combo.index])
  }
  length.na<-sum(is.na(ages.of.first))
  length.all<-length(ages.of.first)
  ages.of.first<-ages.of.first[which(!is.na(ages.of.first))]
   density.result<-density(ages.of.first, from=0, to=136)
  density.result$y<-.3*(1-(length.na/length.all))*density.result$y/max(density.result$y)
  polygon(c(density.result$x, rev(density.result$x)), 9-c(density.result$y+combo.index, combo.index - rev(density.result$y)), col="gray", border=NA)
  text(x=0, y=9-combo.index, labels=combo.index, pos=2)
  text(x=145,y=9-combo.index, labels=round(length.na/length.all,2))
  lines(x=rep(quantile(ages.of.first,0.025),2), y=9-combo.index+c(-.3,.3))
    lines(x=rep(quantile(ages.of.first,.5),2), y=9-combo.index+c(-.3,.3),lwd=2)

  lines(x=rep(quantile(ages.of.first,.975),2), y=9-combo.index+c(-.3,.3))
print(paste(combo.index, quantile(ages.of.first,.975), range(ages.of.first)))
 }
}
