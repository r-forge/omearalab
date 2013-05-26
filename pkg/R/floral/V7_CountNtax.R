library(RColorBrewer)
library(stringr)
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_NewSimulator.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_UtilityFns.R")
dirs<-system("ls -1 | grep old | grep full | grep -v RData | grep -v '\\.R' | grep -iv pdf",intern=TRUE)
#mypalette<-brewer.pal(8,"Dark2")
mypalette<-heat.colors(12,alpha=0.3)
original.data<-read.csv("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv",stringsAsFactors=FALSE)
observed.ntax<-dim(original.data)[1]
all.histories.list<-list()
combo.names<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")

results.df<-c()

observed.proportions<-rep(NA, length(combo.names))
for (i in sequence(length(combo.names))) {
  chars<-strsplit(combo.names[i],"")[[1]]
  small.dataset<-original.data[which(original.data[,1+1]==chars[1]),] #look at the relevant chars
  small.dataset<-small.dataset[which(small.dataset[,3+1]==chars[3]),]
  small.dataset<-small.dataset[which(small.dataset[,4+1]==chars[4]),]
  observed.proportions[i]<-dim(small.dataset)[1]/dim(original.data)[1]
}
names(observed.proportions)<-combo.names

system("rm *.pdf")

for (dir.index in sequence(length(dirs))) {
  ntax.scaling<-as.numeric(str_extract(str_extract("transonly_bd_20000000_rescale_ntax.old.scale_158710", "scale_\\d+"), "\\d+"))
  #print(dirs[dir.index])
  setwd(dirs[dir.index])
  ntax.vector<-c()
#    file.list<-system("ls -1 | grep total151.RSave",intern=TRUE)
  file.list<-system("ls -1 | grep total137.RSave",intern=TRUE)
  if(grepl("old", dirs[dir.index])) {
	file.list<-system("ls -1 | grep RSave", intern=TRUE)
}
  if (length(file.list)>0) {
    file.list<-paste("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013/", dirs[dir.index], "/", file.list, sep="")
  }
  #  if(!grepl("stacey", dirs[dir.index])) {
  #    setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013Additional")
  #    setwd(dirs[dir.index])
  #    file.list2<-system("ls -1 | grep total137.RSave",intern=TRUE)
  #    if (length(file.list2)>0) {
  #      file.list2<-paste("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013Additional/", dirs[dir.index], "/", file.list2, sep="")
  #    }
  #    file.list<-c(file.list, file.list2)
  #  }
  # print(file.list)
  setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013")
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
  results.df<-rbind(results.df, data.frame(ntax.rescale=ntax.rescale, ntax.actual=median(ntax.vector)))
  print(c(ntax.rescale, median(ntax.vector)))
  setwd("..")

}
#system("open All.PDF")
print(rsults.df)
write.csv(results.df, file="~/Dropbox/ntax.results.csv")
