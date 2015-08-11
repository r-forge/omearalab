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
equilibria <- c(0.0024326401, 0.0370846124, 0.0269083723, 0.9005043615, 0.0002350799, 0.0025138161, 0.0014501831, 0.0288709346)

GetFirstAppearance<-function(counts) {
  age<-NA
  nonzeros<-which(counts!=0)
    	#print(paste("length(nonzeros)", length(nonzeros)))

  if(length(nonzeros)>0) {
    age<-nonzeros[1]-1 
  } else {
  	return(NA)	
  }
  #print(page("age",age))
  if(age>136) {
  	age<-NA
  }
  return(age)
}

CountNonzeros<-function(x) {
	return(sum(x>0))
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
plot(x=c(-10,150), y=c(-3,9), yaxt="n", xlab="MY from present", bty="n", type="n",ylab="",xaxt="n")
axis(side=1, at=c(0, 65, 136, 151), labels=c(-136, -65, 0, "+15"))
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
   distances<-c()
  for (i in sequence(length(file.list))) {
    load(file.list[i])
    #print(tail(history))
    #  print(file.list[i])
    proportional.history<-history
    for (row.index in sequence(dim(history)[1])) {
      proportional.history[row.index, 2:9]<-history[row.index, 2:9]/sum( history[row.index, 2:9])
    }
    distances<-append(distances, dist(as.matrix(rbind(observed.proportions, proportional.history[137, 2:9])))[1]/sqrt(2)) #sqrt(2) is the max distance
    #mat.for.dist<<-as.matrix(rbind(observed.proportions, proportional.history[137, 2:9]))
    #dput(as.matrix(rbind(observed.proportions, proportional.history[137, 2:9])))
    #print(dist(as.matrix(rbind(observed.proportions, proportional.history[137, 2:9])))[1])
    #print(dist(as.matrix(rbind(proportional.history[137, 2:9], proportional.history[137, 2:9])))[1])
    #print(rbind(observed.proportions, proportional.history[137, 2:9]))
    if (i==1) {
      combo.names<-colnames(history)[2:9]
      for (j in 2:9) {
        all.histories.list[[j-1]]<-proportional.history[,c(1,j)] 
      }
    }
    else {
      for (j in 2:9) {
        all.histories.list[[j-1]]<-cbind( all.histories.list[[j-1]], proportional.history[,j])
      }
    }
    #print(proportional.history)
    # print(paste("ntax present = ",sum(history[137,2:9])))
    ntax.vector<-append(ntax.vector, sum(history[137,2:9]))
    rm(history)
    rm(proportional.history)
  }

  for(combo.index in rev(sequence(8))) {
  	     current<-all.histories.list[[combo.index]]
		
  		ages.of.first<-c()
  		for (i in sequence(length(file.list))) {
	    	load(file.list[i])
    		ages.of.first<-append(ages.of.first, GetFirstAppearance(history[,combo.index+1]))
  		}
	  length.na<-sum(is.na(ages.of.first))
	  length.all<-length(ages.of.first)
	#  ages.of.first<-ages.of.first[which(!is.na(ages.of.first))]

#	  density.result<-density(ages.of.first, from=0, to=136)
#	  density.result$y<-.3*(1-(length.na/length.all))*density.result$y/max(density.result$y)
#	    lines(x=c(0,136),y=rep(9-combo.index,2), col="gray")
		times<-c(current[,1], rev(current[,1]))
		median.proportions<-apply(current[,-1],1,quantile,probs=0.5)
		ci.0.975.proportions<-apply(current[,-1],1,quantile,probs=0.975)
		ci.0.75.proportions<-apply(current[,-1],1,quantile,probs=0.75)
		ci.0.025.proportions<-apply(current[,-1],1,quantile,probs=0.025)
		ci.0.25.proportions<-apply(current[,-1],1,quantile,probs=0.25)
		ci.0.proportions<-apply(current[,-1],1,quantile,probs=0)
		ci.1.proportions<-apply(current[,-1],1,quantile,probs=1)

		offset<-8.5

		print(length(times))
		print(length(median.proportions))
		print(range(offset-combo.index+c(median.proportions, rev(median.proportions))))
#print(length(offset-c(apply(current[,-1],1,quantile,probs=0.5)+combo.index, combo.index - rev(apply(current[,-1],1,quantile,probs=0.5)))))
		fill.color<-rgb(0,0,0,0.1)
#	  polygon(times, offset-combo.index - 0.5 + c(ci.0.975.proportions, 0*rev(ci.0.975.proportions)), col= fill.color, border="NA")
#	  polygon(times, offset-combo.index - 0.5 + c(ci.0.75.proportions, 0*rev(ci.0.75.proportions)), col= fill.color, border="NA")
#	  polygon(times, offset-combo.index - 0.5 + c(median.proportions, 0*rev(median.proportions)), col=fill.color, border="NA")
	#  polygon(times, offset-combo.index - 0.5 + c(ci.0.025.proportions, 0*rev(ci.0.025.proportions)), col="pink", border="NA")
		polygon(times, offset-1.3*combo.index - 0.5 + c(ci.1.proportions, rev(ci.0.proportions)), col= fill.color, border="NA")

		polygon(times, offset-1.3*combo.index - 0.5 + c(ci.0.975.proportions, rev(ci.0.025.proportions)), col= fill.color, border="NA")
		polygon(times, offset-1.3*combo.index - 0.5 + c(ci.0.75.proportions, rev(ci.0.25.proportions)), col= fill.color, border="NA")
		lines(current[,1],offset-1.3*combo.index - 0.5 + median.proportions, col="red")


	  #lines(current[,1], offset - combo.index + 0.5*median.proportions)
	  #lines(current[,1], offset - combo.index - 0.5*median.proportions)

	  text(x=0, y=offset-1.3*combo.index, labels=gsub("x","",combo.names[combo.index]), pos=2)
	  #lines(x=rep(quantile(ages.of.first,0.025),2), y=offset-combo.index+c(-.3,.3))
	  #  lines(x=rep(quantile(ages.of.first,.5),2), y=offset-combo.index+c(-.3,.3),lwd=2)

	  #lines(x=rep(quantile(ages.of.first,.975),2), y=offset-combo.index+c(-.3,.3))	
	  lines(x=rep(136,2), y=c(offset-1.3*combo.index +0.5, offset - 1.3*combo.index -0.5), col="white")
	  	  lines(range(current[,1], max(current[,1])+3), rep(offset - 1.3*combo.index - 0.5, 2))
	  lines(range(current[,1], max(current[,1])+3), rep(offset - 1.3*combo.index + 0.5, 2))

	  points(136, offset-1.3*combo.index -0.5 + observed.proportions[combo.index],pch=18, col='red')
	  points(136+17, offset-1.3*combo.index - 0.5 + equilibria[combo.index], pch="-", col='blue',cex=2)
 	}
}
