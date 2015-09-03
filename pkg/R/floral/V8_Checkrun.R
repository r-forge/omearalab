#system("curl http://cran.r-project.org/src/contrib/Archive/deSolve/deSolve_1.9.tar.gz > ~/Desktop/deSolve_1.9.tar.gz")
#system("R CMD INSTALL ~/Desktop/deSolve_1.9.tar.gz")
#library(diversitree)
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/")
source("V8_CommandsLocal.R")
library(foreach)
library(doMC)
library(parallel)
registerDoMC(max(detectCores()-1, 1))


make.musse.modified.for.sim <- function(tree, states, k, sampling.f=NULL, strict=FALSE,
                       safe=FALSE) {
  cache <- diversitree:::make.cache.musse(tree, states, k, sampling.f, strict)
  branches <- diversitree:::make.branches.musse(k, list(safe))
  root.p.vector=rep(0,k)
  root.p.vector[1]=1
  ll.musse <- function(pars, condition.surv=TRUE, root=ROOT.GIVEN,
                       root.p=root.p.vector, intermediates=FALSE) {
    if ( length(pars) != k*(k+1) )
      stop(sprintf("Invalid length parameters (expected %d)",
                   k*(k+1)))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    diversitree:::ll.xxsse(pars, cache, initial.conditions.musse, branches,
             condition.surv, root, root.p, intermediates)
  }

  ll <- function(pars, ...) ll.musse(pars, ...)
  class(ll) <- c("musse", "function")
  attr(ll, "k") <- k
  ll
}


#THIS FUNCTION ADDED BY BCO Brian O'Meara TO DEAL WITH FAILURE IN THE CASE OF EXTRA PARAMS
starting.point.musse.extra <- function(tree, k, q.div=5, yule=FALSE,argnames) {
  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  r <- if  ( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]
  q <- r/q.div
  p <- rep(q,length(argnames))
  names(p) <- argnames
  p[grep("lambda",names(p))]<-pars.bd[1]
  p[grep("mu",names(p))]<-pars.bd[2]
  p
}


doSingleModelFromHisseSim<-function(phy, F, T, D, S=6) {
  startTime<-proc.time()
  nAngiosperms=250000
  sampling.f<-rep(2^S/nAngiosperms,2^S)
  results.vector.all<-c()
  states <- sapply(phy$tip.state, comboAsDecimal, S=S)
  names(states) <- names(phy$tip.state)
#  lik <- make.musse.modified.for.sim(tree=phy, states=states, k=2^S, sampling.f=sampling.f)
  lik <- make.musse(tree=phy, states=states, k=2^S, sampling.f=sampling.f, strict=FALSE)

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