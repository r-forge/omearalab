setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
library(diversitree)
source("V6_UtilityFns.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_CommandsLocal.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/AICweightJan2013")
load(file="testing.model.subset.Rsave")


doSingleModelFromSim<-function(phy, F, T, D, S=6) {
  nAngiosperms=250000
  sampling.f<-rep(2^S/nAngiosperms,2^S)
  results.vector.all<-c()
  lik <- make.musse.modifiedWithRootFixedAt1(tree=phy, states=phy$tip.state, k=2^S, sampling.f=sampling.f)
  print("Made musse ")
  assign("extralist",list(),envir = .GlobalEnv) #naughty
  lik.trans <- modify_transitions(lik, type=T, F=F, S=S,extralist=extralist)
  argnames(lik.trans)
  lik.final <- modify_diversification(lik.trans, type=D, F=F, S=S,extralist=extralist)
  argnames(lik.final)
  p <- starting.point.musse(phy, 2^S)
  if (length(extralist)>0) {
    p <- starting.point.musse.extra(phy,2^S,argnames=argnames(lik.final)) #this is done as otherwise won't get right starting vector
  }
  print("Starting fit")
  fit.final <- find.mle(lik.final,p,method="subplex",hessian=TRUE, verbose=TRUE) #the hessian lets us get standard errors
  print("Found fit")
  print(fit.final)
  fit.final.se<-fit.final
  fit.se<-rep(NA,length(fit.final$par))
  try(fit.se<-sqrt(diag(pseudoinverse(-1*fit.final$hessian)))) #just for extra protection
  fit.final.se$par<-fit.se
  names(fit.final.se$par)<-names(fit.final$par)
  try(names(fit.se)<-paste(names(fit.final$par),".se",sep=""))
  final.matrix<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),fit.final$par),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(fit.final$par))))
  rownames(final.matrix)<-paste("FINAL_",rownames(final.matrix),sep="") #to make it easier to grep
  print(formatC(final.matrix,format="f",digits=30,drop0trailing=TRUE))
  print(paste("FINAL_F ",F,sep=""))
  print(paste("FINAL_T ",T,sep=""))
  print(paste("FINAL_D ",D,sep=""))
  print(paste("FINAL_S ",S,sep=""))
  print(paste("FINAL_filename ",filename,sep=""))
  elapsedTime<-(proc.time()-startTime)[3]
  print(paste("elapsedTime = ",elapsedTime))
  final.matrix.ml<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),elapsedTime,coef(fit.final,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","elapsedTime",names(coef(fit.final,full=TRUE,extra=TRUE)))))
  final.matrix.se<-matrix(c(NA,NA,NA,NA,NA,NA,NA,coef(fit.final.se,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","elapsedTime",names(coef(fit.final.se,full=TRUE,extra=TRUE)))))
  final.matrix.all<-cbind(mle=final.matrix.ml,se=final.matrix.se)
  colnames(final.matrix.all)<-c("mle","se")
  return(final.matrix.all)
}

generating.models<-testing.model.subset[c(1, 2, 51),]
nreps<-1
lambda.columns<-which(grepl("lambda", names(generating.models)))[-1]
mu.columns<-which(grepl("mu", names(generating.models)))[-1]
q.columns<-which(grepl("q", names(generating.models)))[-1]
for (rep in sequence(nreps)) {
  #for (model.index in sequence(dim(generating.models)[1])) {
  for (model.index in sequence(1)) {
      
    pars<-unlist(c(generating.models[model.index, lambda.columns], generating.models[model.index, mu.columns], generating.models[model.index, q.columns]))
    pars<-unlist(c(generating.models[model.index, lambda.columns]-generating.models[model.index, mu.columns], 0*generating.models[model.index, mu.columns], generating.models[model.index, q.columns])) #YULE
    
    names(pars)<-NULL
    phy<-NULL
    while(is.null(phy)) {
      phy<-tree.musse(pars=pars, max.taxa=464, x0=1)
    }
    for (testing.model.index in sequence(1)) {
      
    #for (testing.model.index in sequence(dim(testing.model.subset)[1])) {
      final.matrix.all<-doSingleModelFromSim(phy, F=vectorToString(convertFocalLabelToFocalVector(testing.model.subset$focal[testing.model.index], 6, "x")), T=testing.model.subset$T[testing.model.index], D=testing.model.subset$D[testing.model.index])
    }
    
  }
}