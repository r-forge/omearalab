




source("V6_StochasticSSASims_Functions.R")
load("all.results.cleaned.Rsave")
names(highlevel.dataframe)[1:10]
highlevel.dataframe$deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
highlevel.dataframe$AICweight<-exp( -0.5 * highlevel.dataframe$deltaAIC)
highlevel.dataframe$AICweight <- highlevel.dataframe$AICweight/sum(highlevel.dataframe$AICweight)
highlevel.dataframe<-highlevel.dataframe[order(highlevel.dataframe$AICweight, decreasing=TRUE),]

source("V6_UtilityFns.R")
source("V6_NewSimulator.R")

key.focal.vector<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")



focal.dataframe<-highlevel.dataframe[,]
print("dimensions")
focal.dataframe<-focal.dataframe[which(focal.dataframe$T==1), ]
print(dim(focal.dataframe))
focal.dataframe$deltaAIC<-focal.dataframe$AIC-min(focal.dataframe$AIC)
focal.dataframe$AICweight<-exp( -0.5 * focal.dataframe$deltaAIC)
focal.dataframe$AICweight <- focal.dataframe$AICweight/sum(focal.dataframe$AICweight)
q.values<-matrix(nrow=dim(focal.dataframe)[1],ncol=24)
q.names<-c()
for (focal.index.1 in sequence(length(key.focal.vector))) {
  for (focal.index.2 in sequence(length(key.focal.vector))) {
    if (focal.index.1 != focal.index.2) {
      print(c(focal.index.1,focal.index.2))
      if( vectorMismatch(unlist(strsplit(key.focal.vector[focal.index.1],split="")), unlist(strsplit(key.focal.vector[focal.index.2], split="")))==1) {
        q.names<-append(q.names,paste("q", key.focal.vector[focal.index.1], "_", key.focal.vector[focal.index.2], sep=""))
        #print(paste("q", key.focal.vector[focal.index.1], "_", key.focal.vector[focal.index.2], sep=""))
        matching.indices<-getMatchingTransitionIndices(names(focal.dataframe),  key.focal.vector[focal.index.1], key.focal.vector[focal.index.2])
        q.values[,length(q.names)] <- apply(focal.dataframe[,matching.indices], 1, mean)
      }
    }
  }
}


colnames(q.values)<-q.names

lambda.values<-matrix(nrow=dim(focal.dataframe)[1],ncol=8)
lambda.names<-c()
for (focal.index.1 in sequence(length(key.focal.vector))) {
  lambda.names<-append(lambda.names, paste("lambda",key.focal.vector[focal.index.1],sep=""))
  matching.indices <- grepl(gsub("x",".",paste("lambda",key.focal.vector[focal.index.1],sep="")), names(focal.dataframe), perl=TRUE)
  lambda.values[, length(lambda.names)] <- apply(focal.dataframe[,matching.indices], 1, mean)
}
colnames(lambda.values)<-lambda.names


mu.values<-matrix(nrow=dim(focal.dataframe)[1],ncol=8)
mu.names<-c()
for (focal.index.1 in sequence(length(key.focal.vector))) {
  mu.names<-append(mu.names, paste("mu",key.focal.vector[focal.index.1],sep=""))
  matching.indices <- grepl(gsub("x",".",paste("mu",key.focal.vector[focal.index.1],sep="")), names(focal.dataframe), perl=TRUE)
  mu.values[, length(mu.names)] <- apply(focal.dataframe[,matching.indices], 1, mean)
}
colnames(mu.values)<-mu.names










q.means <- apply(q.values, 2, weighted.mean, w=focal.dataframe$AICweight)
lambda.means <- apply(lambda.values, 2, weighted.mean, w=focal.dataframe$AICweight)
mu.means <- apply(mu.values, 2, weighted.mean, w=focal.dataframe$AICweight)

all.values<-c(q.means, lambda.means,mu.means)
print("all rates")
print(all.values)

#go Yule just to try
#yule.scale<-1 #if 0, no different from regular birth death; if 1, pure birth
#lambda.means<-lambda.means-(yule.scale*mu.means)
#mu.means<-mu.means*(1-yule.scale)

x0 <- c(2, 0, 0, 0, 0, 0, 0, 0) #start with two individuals in state 000000, or at least 0x00xx
names(x0) <- key.focal.vector

parms<-c(q.means, lambda.means, mu.means)
names(parms)


a<-sapply(names(q.means), extractFromAndPaste)
a<-append(a, apply(cbind(key.focal.vector,names(lambda.means)), 1, paste, sep="", collapse="*"))
a<-append(a, apply(cbind(key.focal.vector,names(mu.means)), 1, paste, sep="", collapse="*"))
names(a)<-NULL

#states are rows, transitions are columns
nu<-matrix(0,nrow=length(x0),ncol=length(a))
for (i in sequence(length(a))) {
  if(grepl("lambda",a[i])) {
     combo<-substr(a[i],1,6)
     nu[which(names(x0)==combo),i]<-1
  }
  if(grepl("mu",a[i])) {
     combo<-substr(a[i],1,6)
     nu[which(names(x0)==combo),i]<- (-1)
  }
  if(grepl("q",a[i])) {
    nu[which(names(x0)==extractFrom(a[i])),i]<- (-1)
    nu[which(names(x0)==extractTo(a[i])),i]<- (1)
  }
}

#now to get rid of the x, which the ssa fn does not like

save(list=c("x0", "q.means", "lambda.means", "mu.means"), file="Rates.Rsave",compress=TRUE)
