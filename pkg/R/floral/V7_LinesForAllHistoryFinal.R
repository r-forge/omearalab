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

system("rm *.pdf")

subsample.proportion<-function(p,observed.ntax) {
  recovered.p<-rbinom(1,observed.ntax,p)/observed.ntax
  return(recovered.p)
}
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
  setwd("..")
  print(dirs[dir.index])
  print(paste("number of completed runs=", length(file.list)))
  print(round(quantile(ntax.vector)))
  print("distances of observed freq from simulated freq")
  print(round(quantile(distances), 2))
  if (length(file.list)>1) {
    pdf(file=paste(dirs[dir.index],".pdf",sep=""),height=4)
    par(mfcol=c(2,4))
    proportions.mean<-rep(NA,8)
    for (i in sequence(8)) {
      current<-all.histories.list[[i]]
      plot(x=c(min(current[,1]),1.3*max(current[,1])),y=c(0,1),type="n",bty="n",xaxt="n",yaxt="n", xlab="MY from root",ylab="Proportion",main=paste(combo.names[i],"\n",sub("_bd_20000000_rescale", "", dirs[dir.index]),sep=""))
      polygon(x=c(136, 136+15, 136+15, 136), y=c(0, 0, 1, 1), border=NA, col="gray")
      axis(side=1,at=c(0,136),tick=FALSE)
      axis(side=1,at=seq(from=0,to=136,length.out=5),labels=FALSE)
      axis(side=2,at=c(0,1),tick=FALSE)
      axis(side=2,at=seq(from=0,to=1,length.out=5),labels=FALSE)
      for (j in 2:dim(all.histories.list[[i]])[2]) {
        lines(x=current[,1],y=current[,j],col=mypalette[i])
      }
      # print(dim(current))
      # if (dim(current)[2]==9) {
      try(lines(x=current[,1],y=apply(current[,-1],1,quantile,probs=0.975), col="black",lty="dotted"))
      try(lines(x=current[,1],y=apply(current[,-1],1,quantile,probs=0.025), col="black",lty="dotted"))
      try(lines(x=current[,1],y=apply(current[,-1],1,quantile,probs=0.5), col="black"))
      
      densitylines<-density(current[dim(current)[1],2:dim(current)[2]],from=0, to=1)
      densitylines.subsample<-density(sapply(current[dim(current)[1],2:dim(current)[2]], subsample.proportion, observed.ntax=observed.ntax) ,from=0, to=1)
      polygon(x=10+max(all.histories.list[[i]][,1])+c(0,(10*densitylines$y/max(densitylines$y))), y=c(densitylines$x,0),col=mypalette[i],border=mypalette[i])
      polygon(x=30+max(all.histories.list[[i]][,1])+c(0,(10*densitylines.subsample$y/max(densitylines.subsample$y))), y=c(densitylines.subsample$x,0),col=mypalette[i],border=mypalette[i])
      chars<-strsplit(combo.names[i],"")[[1]]
      small.dataset<-original.data[which(original.data[,1+1]==chars[1]),] #look at the relevant chars
      small.dataset<-small.dataset[which(small.dataset[,3+1]==chars[3]),]
      small.dataset<-small.dataset[which(small.dataset[,4+1]==chars[4]),]
      lines(x=10+max(all.histories.list[[i]][,1])+c(0,11), y=rep(dim(small.dataset)[1]/dim(original.data)[1],2),col="black",lwd=2)
      lines(x=30+max(all.histories.list[[i]][,1])+c(0,11), y=rep(dim(small.dataset)[1]/dim(original.data)[1],2),col="black",lwd=2)
      proportions.mean[i]<-mean(current[137, -1])
    }
    names(proportions.mean)<-combo.names
    summary.df<-rbind(summary.df, cbind(data.frame(model=dirs[dir.index], nrun=length(file.list), ntax.025=round(quantile(ntax.vector,.025)), ntax.med=round(quantile(ntax.vector,.5)), ntax.975=round(quantile(ntax.vector,.975)), dist.025=round(quantile(distances, 0.025), 2), dist.med=round(quantile(distances, 0.5), 2), dist.975=round(quantile(distances, 0.975), 2)), data.frame(t(round(proportions.mean, 2)))))
    
    #}
    dev.off()
  }
}
system('/System/Library/Automator/Combine\\ PDF\\ Pages.action/Contents/Resources/join.py -o ~/Desktop/Sims.pdf *.pdf')
system("mv ~/Desktop/Sims.pdf ~/Dropbox/All.PDF")
system("open ~/Dropbox/All.PDF")
row.names(summary.df)<-NULL
print(summary.df)
write.csv(summary.df, file="~/Dropbox/summary.csv")
