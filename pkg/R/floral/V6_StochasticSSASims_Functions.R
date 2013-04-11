

source("V6_UtilityFns.R")
source("V6_NewSimulator.R")


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
