source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_StochasticSSASims_Functions.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_UtilityFns.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_NewSimulator.R")


CreateRatesFile <- function(constraint="full", net.div=FALSE, x0=NULL, x0.rescale=NULL, AIC.cutoff=20, q.rescale=1, best.only=FALSE) {
  load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/all.results.cleaned.Rsave")
  names(highlevel.dataframe)[1:10]
  highlevel.dataframe$deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
  highlevel.dataframe$AICweight<-exp( -0.5 * highlevel.dataframe$deltaAIC)
  highlevel.dataframe$AICweight <- highlevel.dataframe$AICweight/sum(highlevel.dataframe$AICweight)
  highlevel.dataframe<-highlevel.dataframe[order(highlevel.dataframe$AICweight, decreasing=TRUE),]
  
  
  key.focal.vector<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")
  
  
  
  focal.dataframe<-highlevel.dataframe[,]
  print("dimensions original")
  print(dim(focal.dataframe))
  if(constraint=="transonly") {
    focal.dataframe1<-focal.dataframe[which(focal.dataframe$D==1),]
    focal.dataframe3<-focal.dataframe[which(focal.dataframe$D==3),]
    focal.dataframe<-rbind(focal.dataframe1, focal.dataframe3)
  }
  if(constraint=="divonly") {
    focal.dataframe<-focal.dataframe[which(focal.dataframe$T==1), ]
  }
  if(constraint=="symmetry") {
    focal.dataframe0<-focal.dataframe[which(focal.dataframe$focal=="xx0xxx"),]
    focal.dataframe1<-focal.dataframe[which(focal.dataframe$focal=="xx1xxx"),]
    focal.dataframe<-rbind(focal.dataframe0, focal.dataframe1)
  }
  if(constraint=="staceyfull") {
  	focal.dataframe<-focal.dataframe[which(focal.dataframe$D!=1),]
  	focal.dataframe<-focal.dataframe[which(focal.dataframe$D!=3),]
    focal.dataframe<-focal.dataframe[which(focal.dataframe$T!=1), ]
  }
  if(constraint=="transonlyRedo") {
    focal.dataframe1<-focal.dataframe[which(focal.dataframe$D==1),]
    focal.dataframe3<-focal.dataframe[which(focal.dataframe$D==3),]
    focal.dataframe<-rbind(focal.dataframe1, focal.dataframe3)
  }

  
  print("dimensions post filter")
  print(dim(focal.dataframe))
  print("table T")
  print(table(focal.dataframe$TransitionModel))
  print("table D")
  print(table(focal.dataframe$DiversificationModel))
  print("table focal")
  print(table(focal.dataframe$focal))
  
  
  focal.dataframe$deltaAIC<-focal.dataframe$AIC-min(focal.dataframe$AIC)
  print("dim focal.dataframe before AIC cutoff")
  print(dim(focal.dataframe))
  focal.dataframe<-focal.dataframe[which(focal.dataframe$deltaAIC<AIC.cutoff),]
  if(best.only) {
    focal.dataframe<-focal.dataframe[which(focal.dataframe$deltaAIC<=0) ,]
    if(dim(focal.dataframe)[1]!=1) {
      stop("dimension of focal.dataframe not correct for using best model") 
    }
  }
  print("dim focal.dataframe after AIC cutoff")
  print(dim(focal.dataframe))
  focal.dataframe$AICweight<-exp( -0.5 * focal.dataframe$deltaAIC)
  focal.dataframe$AICweight <- focal.dataframe$AICweight/sum(focal.dataframe$AICweight)
  q.values<-matrix(nrow=dim(focal.dataframe)[1],ncol=24)
  q.names<-c()
  for (focal.index.1 in sequence(length(key.focal.vector))) {
    for (focal.index.2 in sequence(length(key.focal.vector))) {
      if (focal.index.1 != focal.index.2) {
        #print(c(focal.index.1,focal.index.2))
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
  
  diversification.values <- lambda.values - mu.values
  diversification.names <- c()
  for (focal.index.1 in sequence(length(key.focal.vector))) {
    diversification.names<-append(diversification.names, paste("div",key.focal.vector[focal.index.1],sep=""))
  }
  colnames(diversification.values)<-diversification.names
  
  
  
  
  
  
  
  
  
  q.means <- q.rescale*apply(q.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  ef.means <- apply(mu.values/lambda.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  turnover.means <- apply(mu.values + lambda.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  netdiv.means <- apply(lambda.values - mu.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  #lambda.means <- apply(lambda.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  #mu.means <- apply(mu.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  lambda.means <- rep(NA, length(ef.means))
  names(lambda.means)<-lambda.names
  mu.means <- rep(NA, length(ef.means))
  names(mu.means)<-mu.names
  for (i in sequence(length(ef.means))) {
    # lambda.means[i] <- getB.ef(ef=ef.means[i], turn=turnover.means[i])
    #  mu.means[i] <- getD.ef(ef=ef.means[i], turn=turnover.means[i])
    lambda.means[i] <- getB.net(net.div=netdiv.means[i], turn=turnover.means[i])
    mu.means[i] <- getD.net(net.div=netdiv.means[i], turn=turnover.means[i])
    
  }
  #diversification.means <- apply(diversification.values, 2, weightedHarmonicMeanZeroCorrection, w=focal.dataframe$AICweight)
  diversification.means <- lambda.means - mu.means
  names(diversification.means) <- diversification.names
  
  all.values<-c(q.means, lambda.means, mu.means, diversification.means)
  print("all rates")
  print(all.values)
  
  if(net.div) {
    lambda.means<-diversification.means
    names(lambda.means)<-lambda.names
    mu.means<-0*mu.means
  }
  
  #go Yule just to try
  #yule.scale<-1 #if 0, no different from regular birth death; if 1, pure birth
  #lambda.means<-lambda.means-(yule.scale*mu.means)
  #mu.means<-mu.means*(1-yule.scale)
  
  if(is.null(x0)) {
    x0 <- c(2, 0, 0, 0, 0, 0, 0, 0) #start with two individuals in state 000000, or at least 0x00xx
    names(x0) <- key.focal.vector
  }
  if(is.null(x0.rescale)) {
    x0.rescale<-x0
  }
  
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
  
  save(list=c("x0", "q.means", "lambda.means", "mu.means", "diversification.means", "nu", "a", "constraint", "net.div", "x0.rescale"), file="Rates.Rsave",compress=TRUE)
}

#constraint can be "full", "transonly", "divonly", "symmetry"
MakeRunFiles<-function(constraint="full", net.div=FALSE, x0=NULL, x0.rescale=NULL, tf=156, t.rescale=136, submit=FALSE, nrep=50, q.rescale=1, best.only=FALSE) {
  file.string<-constraint
  if (net.div) {
    file.string<-paste(file.string, "netdiv", sep="_") 
  } else {
    file.string<-paste(file.string, "bd", sep="_")
  }
  if (!is.null(x0)) {
    file.string<-paste(file.string, vectorToString(x0), sep="_")
  }
  if (!is.null(x0.rescale)) {
    file.string<-paste(file.string, "rescale", sep="_")
  }
  if (best.only) {
    file.string<-paste(file.string, "best", sep="_") 
  }
  if (q.rescale!=1) {
    file.string<-paste(file.string, "q.rescale", q.rescale, sep="_") 
  }
  
  system(paste("mkdir ",file.string))
  setwd(file.string)
  system("cp /Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7*.R .")
  CreateRatesFile(constraint=constraint, net.div=net.div, x0=x0, x0.rescale=x0.rescale, q.rescale=q.rescale, best.only=best.only)
  
  cat('if(!require(gmp)) {
  install.packages("gmp",repos="http://cran.us.r-project.org", lib=tempdir())
  .libPaths(c(.libPaths(), tempdir()))
}
source("V7_UtilityFns.R")
source("V7_NewSimulator.R")
source("V7_StochasticSSASims_Functions.R")
load("Rates.Rsave")
      
', file="ActualRun.R")
  cat(paste("doParallelSSA(tf=",tf,", x0=x0, q.means=q.means, lambda.means=lambda.means, mu.means=mu.means, maxWallTime=Inf, file.string='", file.string, "', rescale.species=250000, yule.scale=0, full.history=FALSE,print.freq=10000, t.rescale=", t.rescale, ", x0.rescale=x0.rescale)", sep=""), file="ActualRun.R", append=TRUE)
  
  cat(paste('#! /bin/sh


echo "I am process id $$ on" `hostname`
/usr/bin/R CMD BATCH --no-save ActualRun.R ', file.string, '.$$.`hostname`.Rout
ls
cat *Rout
pwd
echo "Now finished"
', sep=""), file="myProg")
  system("chmod u+x myProg")
  cat(paste('executable=myProg
universe=vanilla
arguments=', file.string, '$(Cluster).$(Process) 5
requirements = Memory >= 0
output=results.output.$(Process)
error=results.error.$(Process)
transfer_input_files=Rates.Rsave,V7_NewSimulator.R,V7_StochasticSSASims_Functions.R,V7_UtilityFns.R,ActualRun.R
log=results.log
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT_OR_EVICT
queue ', nrep, '

', sep=""), file="myJob.submit")
  if(submit) {
    system("/condor/condor-installed/bin//condor_submit myJob.submit") 
  }
  setwd("..")
}