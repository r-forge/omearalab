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

MakeRunFiles(constraint="full", net.div=FALSE, tf=136, submit=TRUE, nrep=50, best.only=FALSE, x0.rescale=x0.rescale, x0=x0, ntax.old.scale=40348)
MakeRunFiles(constraint="transonly", net.div=FALSE, tf=136, submit=TRUE, nrep=50, best.only=FALSE, x0.rescale=x0.rescale, x0=x0, ntax.old.scale=158710)
MakeRunFiles(constraint="divonly", net.div=FALSE, tf=136, submit=TRUE, nrep=50, best.only=FALSE, x0.rescale=x0.rescale, x0=x0, ntax.old.scale=29042)
ntax.v<-round(exp(seq(from=log(10000), to=log(400000), length.out=25)))
for (i in sequence(length(ntax.v))) {
	MakeRunFiles(constraint="full", net.div=FALSE, tf=136, submit=TRUE, nrep=2, best.only=FALSE, x0.rescale=x0.rescale, x0=x0, ntax.old.scale=ntax.v[i])

}
