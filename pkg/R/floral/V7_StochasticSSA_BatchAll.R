setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013")
system("cp /Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7*R .")
source("V7_StochasticSSASims_CreateAllFilesForRun.R")
library(compiler)
enableJIT(3)
original.data<-read.csv("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv",stringsAsFactors=FALSE)

combo.names<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")

#use tip state freq to figure out 
x0.rescale<-rep(0, 8)
for (i in sequence(8)) {
    chars<-strsplit(combo.names[i],"")[[1]]
    small.dataset<-original.data[which(original.data[,1+1]==chars[1]),] #look at the relevant chars
    small.dataset<-small.dataset[which(small.dataset[,3+1]==chars[3]),]
    small.dataset<-small.dataset[which(small.dataset[,4+1]==chars[4]),]
    x0.rescale[i]<-dim(small.dataset)[1]/dim(original.data)[1]
}

x0<-c(2, 0, 0, 0, 0, 0, 0, 0)
names(x0)<-combo.names
names(x0.rescale)<-combo.names

constraint.vector <- c("full", "transonly", "divonly", "staceyfull")
net.div.vector <- c(FALSE)
best.vector <- c(TRUE, FALSE)
for (i in sequence(length(net.div.vector))) {
	for (j in sequence(length(constraint.vector))) {
    for (k in sequence(length(best.vector)))
		MakeRunFiles(constraint=constraint.vector[j], net.div=net.div.vector[i], tf=65, submit=TRUE, nrep=25, best.only=best.vector[k], x0.rescale=x0.rescale, x0=x0)
	}
}

MakeRunFiles(constraint="full", net.div=FALSE, tf=65, submit=TRUE, nrep=25, best.only=TRUE, x0.rescale=x0.rescale, x0=x0, q.rescale=10)
MakeRunFiles(constraint="full", net.div=FALSE, tf=65, submit=TRUE, nrep=25, best.only=FALSE, x0.rescale=x0.rescale, x0=x0, q.rescale=10)
