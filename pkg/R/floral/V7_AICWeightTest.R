setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
library(diversitree)
source("V6_UtilityFns.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_CommandsLocal.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/AICweightJan2013")
load(file="testing.model.subset.Rsave")
library(foreach)
library(doMC)
registerDoMC(13)

doSingleModelFromSim<-function(phy, F, T, D, S=6) {
  startTime<-proc.time()
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
  fit.final <- find.mle(lik.final,p,method="subplex",hessian=TRUE) #the hessian lets us get standard errors
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
  elapsedTime<-(proc.time()-startTime)[3]
  print(paste("elapsedTime = ",elapsedTime))
  final.matrix.ml<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),elapsedTime,coef(fit.final,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","elapsedTime",names(coef(fit.final,full=TRUE,extra=TRUE)))))
  final.matrix.se<-matrix(c(NA,NA,NA,NA,NA,NA,NA,coef(fit.final.se,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","elapsedTime",names(coef(fit.final.se,full=TRUE,extra=TRUE)))))
  final.matrix.all<-cbind(mle=final.matrix.ml,se=final.matrix.se)
  colnames(final.matrix.all)<-c("mle","se")
  return(final.matrix.all)
}

generating.models<-testing.model.subset[c(1, 2, 51, 51),]
generating.models[dim(generating.models)[1], sapply(generating.models[1,], is.numeric)]<-colMeans(testing.model.subset[, sapply(generating.models[1,], is.numeric)])
generating.models[dim(generating.models)[1], c(1, 2, 5, 7, 8, 9)]<-rep(NA, 6)
generating.models[dim(generating.models)[1], 3]<-"average"
nreps<-100
doRunWithinForeach<-function(rep, generating.models, testing.model.subset) {
  lambda.columns<-which(grepl("lambda", names(generating.models)))[-1]
  mu.columns<-which(grepl("mu", names(generating.models)))[-1]
  q.columns<-which(grepl("q", names(generating.models)))[-1]
  
  for (model.index in sequence(dim(generating.models)[1])) {
  #for (model.index in sequence(1)) {
      
    pars<-unlist(c(generating.models[model.index, lambda.columns], generating.models[model.index, mu.columns], generating.models[model.index, q.columns])) #actual
    #pars<-unlist(c(generating.models[model.index, lambda.columns]-generating.models[model.index, mu.columns], 0*generating.models[model.index, mu.columns], generating.models[model.index, q.columns])) #YULE
    #pars<-unlist(c(generating.models[model.index, lambda.columns]-0.999*generating.models[model.index, mu.columns], 0.001*generating.models[model.index, mu.columns], generating.models[model.index, q.columns])) #ALMOST YULE
    
    
    names(pars)<-NULL
    phy<-NULL
    while(is.null(phy)) {
      try(phy<-tree.musse(pars=pars, max.taxa=464, x0=1))
    }
    #for (testing.model.index in sequence(1)) {
      
    for (testing.model.index in sequence(dim(testing.model.subset)[1])) {
      final.matrix.all<-NULL
      try(final.matrix.all<-doSingleModelFromSim(phy, F=vectorToString(convertFocalLabelToFocalVector(testing.model.subset$focal[testing.model.index], 6, "x")), T=testing.model.subset$T[testing.model.index], D=testing.model.subset$D[testing.model.index]))
      if (!is.null(final.matrix.all)) {
       qIndices<-grep("^q\\d",row.names(final.matrix.all),perl=TRUE)
       lambdaIndices<-grep("^lambda\\d",row.names(final.matrix.all),perl=TRUE)
        muIndices<-grep("^mu\\d",row.names(final.matrix.all),perl=TRUE)
        T<-testing.model.subset$T[testing.model.index]
        D<-testing.model.subset$D[testing.model.index]
        tmp.dataframe<-data.frame(testing.model.subset$focal[testing.model.index],T,transitionModels[T,4],D,diversificationModels[D,5],final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),1],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),1]) 
        names(tmp.dataframe)<-c("focal","T","TransitionModel","D","DiversificationModel","lnLik","AIC","k_all","k_q","k_lambda","k_mu")	
        tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(final.matrix.all[qIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[qIndices,1])))),data.frame(matrix(final.matrix.all[lambdaIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[lambdaIndices,1])))),data.frame(matrix(final.matrix.all[muIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[muIndices,1])))))
        save(phy, rep, pars, final.matrix.all, tmp.dataframe, testing.model.index, model.index, file=paste("R",rep,"Generating",model.index,"Used",testing.model.index,"_",format(Sys.time(), "%H.%M.%S.%b%d.%Y"),".Rsave", sep=""), compress=TRUE)
      }
    }
    
  }
}

foreach  (rep=sequence(nreps)) %dopar% doRunWithinForeach(rep, generating.models, testing.model.subset)

