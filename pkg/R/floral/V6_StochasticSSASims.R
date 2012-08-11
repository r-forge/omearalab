library(GillespieSSA)

load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/Highlevel.dataframe.withrates.Rsave")
names(highlevel.dataframe)[1:10]
highlevel.dataframe$deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
highlevel.dataframe$AICweight<-exp( -0.5 * highlevel.dataframe$deltaAIC)
highlevel.dataframe$AICweight <- highlevel.dataframe$AICweight/sum(highlevel.dataframe$AICweight)
highlevel.dataframe<-highlevel.dataframe[order(highlevel.dataframe$AICweight, decreasing=TRUE),]

source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V6_UtilityFns.R")

key.focal.vector<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")

getMatchingTransitions<-function(focalRow, key.focal.from, key.focal.to, change.positions=c(2, 5, 6), average=TRUE) {
  key.focal.from.split<-strsplit(key.focal.from, "")[[1]]
  key.focal.to.split<-strsplit(key.focal.to, "")[[1]]
  matching.values<-c()
  for (i in c(0:1)) {
    key.focal.from.split[change.positions[1]]<-i
    key.focal.to.split[change.positions[1]]<-i   
    for (j in c(0:1)) {
      key.focal.from.split[change.positions[2]]<-j
      key.focal.to.split[change.positions[2]]<-j
      for (k in c(0:1)) {
        key.focal.from.split[change.positions[3]]<-k
        key.focal.to.split[change.positions[3]]<-k
        #print(names(focalRow)[1:25])
        #print(paste("q", paste(key.focal.from.split, sep="", collapse=""), "_", paste(key.focal.to.split, sep="", collapse=""), sep=""))
        matching.values<-append(matching.values, 
              focalRow[ which(names(focalRow)==paste("q", paste(key.focal.from.split, sep="", collapse=""), "_", paste(key.focal.to.split, sep="", collapse=""), sep=""))])
      }
    }
  }
  if (!average) {
    return(matching.values)
  }
  else {
    return(mean(unlist(matching.values)))
  }
}

getMatchingTransitionIndices<-function(name.vector, key.focal.from, key.focal.to, change.positions=c(2, 5, 6)) {
  key.focal.from.split<-strsplit(key.focal.from, "")[[1]]
  key.focal.to.split<-strsplit(key.focal.to, "")[[1]]
  matching.indices<-c()
  for (i in c(0:1)) {
    key.focal.from.split[change.positions[1]]<-i
    key.focal.to.split[change.positions[1]]<-i   
    for (j in c(0:1)) {
      key.focal.from.split[change.positions[2]]<-j
      key.focal.to.split[change.positions[2]]<-j
      for (k in c(0:1)) {
        key.focal.from.split[change.positions[3]]<-k
        key.focal.to.split[change.positions[3]]<-k
        #print(names(focalRow)[1:25])
        #print(paste("q", paste(key.focal.from.split, sep="", collapse=""), "_", paste(key.focal.to.split, sep="", collapse=""), sep=""))
        matching.indices<-append(matching.indices,which(name.vector==paste("q", paste(key.focal.from.split, sep="", collapse=""), "_", paste(key.focal.to.split, sep="", collapse=""), sep="")))
      }
    }
  }
  return(matching.indices)
  
}

focal.dataframe<-highlevel.dataframe[,]

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

#go Yule just to try
#yule.scale<-0 #if 0, no different from regular birth death; if 1, pure birth
#lambda.means<-lambda.means-(yule.scale*mu.means)
#mu.means<-mu.means*(1-yule.scale)

x0 <- c(2, 0, 0, 0, 0, 0, 0, 0) #start with two individuals in state 000000, or at least 0x00xx
names(x0) <- key.focal.vector

parms<-c(q.means, lambda.means, mu.means)
names(parms)

extractFromAndPaste<-function(transition.name) {
  m <- regexpr("q(.+)\\_",transition.name, perl=F)
   return(paste(sub("_","",sub("q","",regmatches(transition.name, m))), transition.name, sep="*"))
}

extractFrom<-function(transition.name) {
  m <- regexpr("q(.+)\\_",transition.name, perl=F)
  return(sub("_","",sub("q","",regmatches(transition.name, m))))
}

extractTo<-function(transition.name) {
  m <- regexpr("\\_(.+)",transition.name, perl=F)
  return(sub("_","",regmatches(transition.name, m)))
}

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



trySim<-function(x0, a, nu, parms, tf=130, maxWallTime=10) {
  names(parms)<-gsub("x","",names(parms))
  names(x0)<-gsub("x","",names(x0))
  a<-gsub("x","",a)
  min.val<-(-1)
  out<-c()
  try.count<-0
while (min.val<0) {
  try.count<-try.count+1
  print(try.count)
  out<-ssa(x0,a,nu,parms,tf=tf,maxWallTime=maxWallTime,method="D",ignoreNegativeState=TRUE)
  min.val<-min(out$data)
    ssa.plot(out)
  print(c(max(out$data[,1]),max(apply(out$data[,2:9],1, sum))))

  #plot(out$data[,1],apply(out$data[,2:9],1, sum))
  }
  print(tail(out$data))
  ssa.plot(out)
  return(out)
}

#out<-trySim(x0, a, nu, parms, tf=130)

probExtinction<-function(lambda, mu, t, N0=2) {
  #after magallon and sanderson
  epsilon<-mu/lambda
  r<-lambda-mu
  beta<-(exp(r*t) - 1) / (exp(r*t) - epsilon)
  alpha<-epsilon * beta
  return(alpha^N0)
}
