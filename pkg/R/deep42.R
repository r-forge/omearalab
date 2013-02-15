##DEEP42 -- Diversification Estimates Extracted from Phylogenies using 42 (err 90) models##

#Overall goal is to have turnover rate or extinction fraction vary as a result of various
#   parameters.
#Note that with logistic growth, it makes the most sense to pass in split.times as actual
#   split.times (including those with no descendants, but in practice, split.times reconstructed
#   on the phylogeny will probably be used instead. Rabosky and Lovette 2008, in a response,
#   argue this is fairly okay (the two measures are correlated).

library(picante)
library(ape)
library(nloptr)
library(parallel)
source("deepSim.R")

Deep<-function(phy, f=1, model=c("yule", "bd"), turnover.logistic=FALSE, eps.logistic=FALSE, turnover.exp=FALSE, eps.exp=FALSE, turnover.inherit=FALSE, eps.inherit=FALSE, turnover.slice=FALSE, eps.slice=FALSE, turnover.ratch=FALSE, eps.ratch=FALSE, n.cores=NULL, nsims=1000, sanity.cutoff=Inf) {
	
	#Helps keep the numbers consistent so we can more confidently assess which interval we are in:	
	options(digits=10)
	#
	obj <- NULL
	#Sets the main parameters to be used in the model:
	if (is.character(model)) {
		if(model=="yule"){
			turnover.param = TRUE
			eps.param = FALSE
			eps.logistic=FALSE
			eps.inherit=FALSE
			eps.slice=FALSE
			eps.exp=FALSE
			eps.ratch=FALSE
		}
		if(model=="bd"){
			turnover.param = TRUE
			eps.param = TRUE
		}
	}
	#There is a method to this madness -- eventually will allow interactive fixing of various parameters while still estimating others:
	if(turnover.inherit==TRUE | turnover.ratch==TRUE){
		turnover.inherit1=TRUE
		turnover.inherit2=TRUE
		turnover.inherit3=TRUE
		turnover.inherit4=TRUE
	}
	if(turnover.inherit==FALSE & turnover.ratch==FALSE){
		turnover.inherit1=FALSE
		turnover.inherit2=FALSE
		turnover.inherit3=FALSE
		turnover.inherit4=FALSE
	}
	if(eps.inherit==TRUE | eps.ratch==TRUE){
		eps.inherit1=TRUE
		eps.inherit2=TRUE
		eps.inherit3=TRUE
		eps.inherit4=TRUE
	}
	if(eps.inherit==FALSE & eps.ratch==FALSE){
		eps.inherit1=FALSE
		eps.inherit2=FALSE
		eps.inherit3=FALSE
		eps.inherit4=FALSE
	}	
	#Sets two important inputs:
	split.times=sort(branching.times(phy), decreasing=TRUE)
	tot_time=max(branching.times(phy))
	
	#Makes a vector of parameters that are going to be estimated:
	pars=c(turnover.param, turnover.inherit1, eps.param, eps.inherit1, turnover.logistic, eps.logistic, turnover.inherit2, turnover.inherit3, turnover.inherit4, eps.inherit2, eps.inherit3, eps.inherit4, turnover.slice, eps.slice, turnover.exp, eps.exp, turnover.ratch, turnover.ratch, eps.ratch, eps.ratch)
	#An index for use later:
	tmp<-pars
	#The trues specify the number of parameters:
	pars[pars==T] <- 1:length(pars[pars==T])
	np<-max(pars)
	#All parameters not estimated are set to 1+max number of estimated parameters:
	pars[pars==0]<-max(pars)+1
	#Function used for optimizing parameters:
	DevOptimize <- function(p, pars, phy, tot_time, f, turnover.inherit, turnover.weight.logistic, turnover.exp, turnover.ratch, eps.inherit, eps.weight.logistic, eps.exp, eps.ratch, split.times, quantile.set, n.cores) {
		
		#Generates the final vector with the appropriate parameter estimates in the right place:
		model.vec <- numeric(length(pars))
		model.vec[] <- c(p, 0)[pars]
		#Now set entries in model.vec to the defaults if we are not estimating them, or, if a vector of values are supplied, make these fixed values
		#in the model. Will be useful for doing likelihood profiles, contours, etc.:
		if(turnover.inherit == FALSE) {
			model.vec[2]=model.vec[1]
		}
		if(eps.inherit == FALSE) {
			model.vec[4]=model.vec[3]
		}
		if(turnover.weight.logistic == TRUE) {
			turnover.weight.logistic=1
		}
		if(turnover.weight.logistic == FALSE) {
			turnover.weight.logistic=0
		}				
		if(eps.weight.logistic == TRUE) {
			eps.weight.logistic=1
		}
		if(eps.weight.logistic == FALSE) {
			eps.weight.logistic=0
		}		
		if(model.vec[5] == 0) {
			model.vec[5] = 1
		}
		if(model.vec[6] == 0) {
			model.vec[6] = 1
		}
		if(pars[5]<(length(p)+1) & abs(model.vec[5])<Ntip(phy) | pars[6]<(length(p)+1) & abs(model.vec[6])<Ntip(phy)){
				return(10000000)
		}
		#An input of PhiAll:
		turnover.splits <- exp(rnorm(length(split.times), log(model.vec[2]), model.vec[13]))
		eps.splits <- exp(rnorm(length(split.times), log(model.vec[4]), model.vec[14]))
		#
		if(pars[9] < max(pars) | pars[13] < max(pars) | pars[12] < max(pars) | pars[14] < max(pars)){
			if(is.null(n.cores)){
				tot.logl <-c()
				for(i in 1:length(quantile.set)){
					tot.logl <- c(tot.logl,GetLikelihood(phylo=phy, tot_time, f=f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[13], turnover.weight.anc.0=model.vec[7], turnover.weight.anc.half=model.vec[8], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[15], turnover.sigma.anc=model.vec[9], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[14], eps.weight.anc.0=model.vec[10], eps.weight.anc.half=model.vec[11], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[16], eps.sigma.anc=model.vec[12], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], turnover.prob.kick=model.vec[17], turnover.kick.value=model.vec[18], eps.prob.kick=model.vec[19], eps.kick.value=model.vec[20], scaled.set=quantile.set[[i]]))
				}
				phi<-PhiAll(max.time=max(branching.times(phy)), f=f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[13], turnover.weight.anc.0=model.vec[7], turnover.weight.anc.half=model.vec[8], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[15], turnover.sigma.anc=model.vec[9], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[14], eps.weight.anc.0=model.vec[10], eps.weight.anc.half=model.vec[11], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[16], eps.sigma.anc=model.vec[12], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], turnover.prob.kick=model.vec[17], turnover.kick.value=model.vec[18], eps.prob.kick=model.vec[19], eps.kick.value=model.vec[20], turnover.splits, eps.splits, t.edge=0, precision, sanity.cutoff=Inf, n.cores=n.cores)
				scaled.logl<-tot.logl[is.finite(tot.logl)] - log(1-phi)
				logl<-mean(scaled.logl,na.rm=TRUE)
				if(!is.finite(logl)){
					return(10000000)
				}
			}
			else{
				SimSigma<-function(nstart){
					tmp = matrix(,1,ncol=1)
					tot.logl <- GetLikelihood(phylo=phy, tot_time, f=f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[13], turnover.weight.anc.0=model.vec[7], turnover.weight.anc.half=model.vec[8], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[15], turnover.sigma.anc=model.vec[9], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[14], eps.weight.anc.0=model.vec[10], eps.weight.anc.half=model.vec[11], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[16], eps.sigma.anc=model.vec[12], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], turnover.prob.kick=model.vec[17], turnover.kick.value=model.vec[18], eps.prob.kick=model.vec[19], eps.kick.value=model.vec[20], scaled.set=quantile.set[[nstart]])
					tmp[,1] = tot.logl
					tmp
				}
				sim.set<-mclapply(1:length(quantile.set), SimSigma, mc.cores=n.cores)
				tot.logl<-unlist(sim.set)
				phi<-PhiAll(max.time=max(branching.times(phy)), f=f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[13], turnover.weight.anc.0=model.vec[7], turnover.weight.anc.half=model.vec[8], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[15], turnover.sigma.anc=model.vec[9], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[14], eps.weight.anc.0=model.vec[10], eps.weight.anc.half=model.vec[11], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[16], eps.sigma.anc=model.vec[12], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], turnover.prob.kick=model.vec[17], turnover.kick.value=model.vec[18], eps.prob.kick=model.vec[19], eps.kick.value=model.vec[20], turnover.splits, eps.splits, t.edge=0, precision, sanity.cutoff=Inf, n.cores=n.cores)
				scaled.logl<-tot.logl[is.finite(tot.logl)] - log(1-phi)
				logl<-mean(scaled.logl,na.rm=TRUE)
			}
			return(-logl)
		}
		else{
			raw.logl <- GetLikelihood(phylo=phy, tot_time, f=f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[13], turnover.weight.anc.0=model.vec[7], turnover.weight.anc.half=model.vec[8], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[15], turnover.sigma.anc=model.vec[9], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[14], eps.weight.anc.0=model.vec[10], eps.weight.anc.half=model.vec[11], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[16], eps.sigma.anc=model.vec[12], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], turnover.prob.kick=model.vec[17], turnover.kick.value=model.vec[18], eps.prob.kick=model.vec[19], eps.kick.value=model.vec[20], scaled.set=quantile.set)
			phi<-PhiAll(max.time=max(branching.times(phy)), f=f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[13], turnover.weight.anc.0=model.vec[7], turnover.weight.anc.half=model.vec[8], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[15], turnover.sigma.anc=model.vec[9], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[14], eps.weight.anc.0=model.vec[10], eps.weight.anc.half=model.vec[11], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[16], eps.sigma.anc=model.vec[12], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], turnover.prob.kick=model.vec[17], turnover.kick.value=model.vec[18], eps.prob.kick=model.vec[19], eps.kick.value=model.vec[20], turnover.splits, eps.splits, t.edge=0, precision, sanity.cutoff=Inf, n.cores=n.cores)
			logl <- raw.logl - log(1-phi)
			if(!is.finite(logl)){
				return(10000000)
			}
			return(-logl)
		}
	}

	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)

	cat("Initializing...", "\n")
	
	#Begin with a rough search to find estimates from constant bd to be used as starting points in thorough search:
	init.quantile.set<-GetQuantiles(phylo=phy)
	init.pars=c(turnover.param, FALSE, eps.param, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
	init.tmp<-init.pars
	init.pars[init.pars==TRUE] <- 1:length(init.pars[init.pars==TRUE])
	init.pars[init.pars==0]<-max(init.pars)+1
	init.set.pars <- c(1,1,0.5,0.5,Ntip(phy)*2,Ntip(phy)*2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0,0,0,0,0)
	init.ip <- init.set.pars[init.tmp==TRUE]
	init.set.lower <- c(0,0,0,0,-1e6,-1e6,0,0,0,0,0,0,0,0,-10,-10,0,0,0,0)
	lower.init <- init.set.lower[init.tmp==TRUE]
	init.set.upper <- c(10,10,10,10,1e6,1e6,1,1,10,1,1,10,10,10,10,10,1,10,1,10)
	upper.init <- init.set.upper[init.tmp==TRUE]
	init = nloptr(x0=init.ip, eval_f=DevOptimize, lb=lower.init, ub=upper.init, opts=opts, pars=init.pars, phy=phy, tot_time=tot_time, f=f, turnover.inherit=FALSE, turnover.weight.logistic=FALSE, turnover.exp=FALSE, turnover.ratch=FALSE,eps.inherit=FALSE, eps.weight.logistic=FALSE, eps.exp=FALSE, eps.ratch=FALSE, split.times=split.times, quantile.set=init.quantile.set, n.cores=n.cores)

	def.set.pars <- c(init$solution[1],init$solution[1],init$solution[2],init$solution[2],Ntip(phy)*2,Ntip(phy)*2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0,0,0,0,0)
	#Set initials using estimates from constant bd model:
	ip <- def.set.pars[tmp==TRUE]
	def.set.lower <- c(0,0,0,0,-1e6,-1e6,0,0,0,0,0,0,0,0,-10,-10,0,0,0,0)
	lower <- def.set.lower[tmp==TRUE]
	def.set.upper <- c(10,10,10,10,1e6,1e6,1,1,10,1,1,10,10,10,10,10,1,10,1,10)
	upper <- def.set.upper[tmp==TRUE]

	#If inheritance model is chosen, create a list of quantiles for each edge in the tree and perform 3 independent analyses:
	if(turnover.inherit==TRUE | eps.inherit==TRUE | turnover.ratch==TRUE | eps.ratch==TRUE){ 
		cat("Finished. Begin thorough search...", "\n")
		GetQuantileSet<-function(nstarts){
			GetQuantiles(phylo=phy)
		}
		quantile.set<-lapply(1:nsims, GetQuantileSet)
		res<-c()
		for(i in 1:3){
		cat("Restart", i, "of 3", "\n")
			out = nloptr(x0=ip, eval_f=DevOptimize, lb=lower, ub=upper, opts=opts, pars=pars, phy=phy, tot_time=tot_time, f=f, turnover.inherit=turnover.inherit, turnover.weight.logistic=turnover.logistic, turnover.exp=turnover.exp, turnover.ratch=turnover.ratch, eps.inherit=eps.inherit, eps.weight.logistic=eps.logistic, eps.exp=eps.exp, eps.ratch=eps.ratch, split.times=split.times, quantile.set=quantile.set, n.cores=n.cores)
			res<-rbind(res,c(-out$objective,out$solution))
		}
		loglik <- res[which.max(res[,1]),1]
		solution <- numeric(length(pars))
		solution[] <- c(res[which.max(res[,1]),2:dim(res)[2]],0)[pars]
	}
	#If no inheritance, perform thorough search:
	else{
		cat("Finished. Begin thorough search...", "\n")
		quantile.set<-GetQuantiles(phylo=phy)
		out = nloptr(x0=ip, eval_f=DevOptimize, lb=lower, ub=upper, opts=opts, pars=pars, phy=phy, tot_time=tot_time, f=f, turnover.inherit=turnover.inherit, turnover.weight.logistic=turnover.logistic, turnover.exp=turnover.exp, turnover.ratch=turnover.ratch, eps.inherit=eps.inherit, eps.weight.logistic=eps.logistic, eps.exp=eps.exp, eps.ratch=eps.ratch, split.times=split.times, quantile.set=quantile.set, n.cores=n.cores)
		loglik <- -out$objective
		#Recreate the model vector for use in the print function:
		solution <- numeric(length(pars))
		solution[] <- c(out$solution, 0)[pars]
	}
	
	obj = list(loglik = loglik, AIC = -2*loglik+2*np,AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par = pars, opts=opts, f=f, phy=phy, iterations=out$iterations) 
	
	class(obj)<-"diversity"		
	return(obj)	
}

SetParameter <- function(stop.time, param.anc, sigma.indep, weight.anc.0, weight.anc.half, weight.logistic, trend.exponent, split.times, k, param.splits, t.edge) {
	position <- max(which(split.times>=stop.time),1)
	n.taxa <- 1+length(which(split.times>=stop.time))
	param.indep <- param.splits[position]
	weight.anc <- weight.anc.0 * exp(-t.edge * log(2) / weight.anc.half) #so if weight.anc.halflife = INF, is the same as our original weight.anc model
	weight.anc[is.na(weight.anc)==TRUE]=0
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc) #same as autoregressive model
	logistic.scaling <- 1 - (weight.logistic * (n.taxa / k))
	param.mean <- (param.starting * (n.taxa ^ trend.exponent)) * logistic.scaling
	return(param.mean)
}

SetBirth <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	turnover <- SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc.0=turnover.weight.anc.0, weight.anc.half=turnover.weight.anc.half, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits, t.edge=t.edge)
	extinction.fraction <- SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc.0=eps.weight.anc.0, weight.anc.half=eps.weight.anc.half, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits, t.edge=t.edge)
	return(SetBirthGivenParameters(turnover, extinction.fraction))
}

SetBirthGivenParameters <- function(turnover, extinction.fraction) {
	birth.expected <- turnover / (1 + extinction.fraction)
	return(birth.expected)
}

SetDeath <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	turnover<-SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc.0=turnover.weight.anc.0, weight.anc.half=turnover.weight.anc.half,weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits, t.edge=t.edge)
	extinction.fraction<-SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc.0=eps.weight.anc.0, weight.anc.half=eps.weight.anc.half, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits, t.edge=t.edge)
	return(SetDeathGivenParameters(turnover, extinction.fraction))
}

SetDeathGivenParameters <- function(turnover, extinction.fraction) {
	death.expected<-(extinction.fraction * turnover) / (1 + extinction.fraction)
	return(death.expected)
}

SetDiversification <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	turnover<-SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc.0=turnover.weight.anc.0, weight.anc.half=turnover.weight.anc.half, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits, t.edge=t.edge)
	extinction.fraction<-SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc.0=eps.weight.anc.0, weight.anc.half=eps.weight.anc.half, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits, t.edge=t.edge)
	return(SetBirthGivenParameters(turnover, extinction.fraction) - SetDeathGivenParameters(turnover, extinction.fraction))
}

######################################################################################################################################
######################################################################################################################################
### This is our solution for integating diversification  
######################################################################################################################################
######################################################################################################################################

IntegrateSetDiversificationExpected <- function(lower, upper, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic,turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {	
	if (turnover.weight.logistic>0 | eps.weight.logistic>0 | turnover.trend.exponent>0 | eps.trend.exponent>0) {
		boundaries<-split.times[which(split.times>lower)]
		boundaries<-boundaries[which(boundaries<upper)]
		lower.bounds<-sort(c(lower, boundaries))
		upper.bounds<-sort(c(boundaries, upper))
		integration.result<-0
		for (i in sequence(length(lower.bounds))) {
			#piecewise equation: (exp(b-d)*t-exp(b-d)*s)*b / (b-d) 
			integration.result <- integration.result + ((exp(SetDiversification(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*upper.bounds[i])-exp(SetDiversification(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*lower.bounds[i]))*SetBirth(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge))/SetDiversification(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)
		}
		return(integration.result)
	}
	else {
		#piecewise equation: (exp(b-d)*t-exp(b-d)*s)*b / (b-d)
		integration.result<-((exp(SetDiversification(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*upper)-exp(SetDiversification(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*lower))*SetBirth(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)) / SetDiversification(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)
	}
	return(integration.result)
}

IntegrateDiversificationOverTime <- function(stop.time, start.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	return(IntegrateSetDiversificationExpected(lower=stop.time, upper=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge))
}

IntegrateDiversificationOverTime.int.int <- function(stop.time, start.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half,turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	result<-IntegrateDiversificationOverTime(stop.time=stop.time, start.time=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)
	return(result)
}

######################################################################################################################################
######################################################################################################################################
### This function computes Phi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

PhiMorlon <- function(stop.time, f, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	morlon.phi <- (1-exp((SetDiversification(stop.time=stop.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*stop.time)-(SetDiversification(stop.time=0, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*0)) / (1/f+IntegrateDiversificationOverTime.int.int(stop.time=0,start.time=stop.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)))
	return(morlon.phi)	
}

PhiAll <- function(max.time, f, turnover.param.root, turnover.param.indep, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, turnover.sigma.anc, eps.param.root, eps.param.indep, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, eps.sigma.anc, split.times, turn.k, eps.k, turnover.prob.kick, turnover.kick.value, eps.prob.kick, eps.kick.value, turnover.splits, eps.splits, t.edge, precision, sanity.cutoff, n.cores) {
	if (turnover.sigma.anc==0 && turnover.weight.anc.0==0 && turnover.weight.anc.half==0 && eps.sigma.anc==0 && eps.weight.anc.0==0 && eps.weight.anc.half==0){
		return(PhiMorlon(stop.time=max.time, f=f, turnover.param.anc=turnover.param.root, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.root, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge))
	}else{
		return(PhiSim(max.time=max.time, f=f, turnover.param.anc=turnover.param.root, turnover.param.indep=turnover.param.indep, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, turnover.sigma.anc=turnover.sigma.anc, eps.param.anc=eps.param.root, eps.param.indep=eps.param.indep, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.sigma.anc=eps.sigma.anc, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.prob.kick=turnover.prob.kick, turnover.kick.value=turnover.kick.value, eps.prob.kick=eps.prob.kick, eps.kick.value=eps.kick.value, turnover.splits=turnover.splits, eps.splits=eps.splits, precision=precision, sanity.cutoff=sanity.cutoff, n.cores=n.cores))
	}
}

PhiSim <- function(max.time, f, turnover.param.anc, turnover.param.indep, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, turnover.sigma.anc, eps.param.anc, eps.param.indep, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, eps.sigma.anc, split.times, turn.k, eps.k, turnover.prob.kick, turnover.kick.value, eps.prob.kick, eps.kick.value, turnover.splits, eps.splits, precision, sanity.cutoff, n.cores) {

	phi.morlon.estimate<-PhiMorlon(stop.time=max.time, f, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=0)
	if (is.na(phi.morlon.estimate)) {
		phi.morlon.estimate<-1-.Machine$double.eps^.5
	}
	#Binomial variance is np(1 − p). So 
	#v = np(1-p)
	#v/((p(1-p))) = n
	#Need to futz:
	precision<-100
	n.sims<-ceiling(precision/(phi.morlon.estimate * (1 - phi.morlon.estimate)))
	#Actually saving seed in order to use a different seed for the 
	save.seed<-round(runif(1,0,100000000))
	set.seed(42)
	if(is.null(n.cores)){
		actual.phi<-sum(unlist(sapply(replicate(n.sims, DeepSim(max.time=max.time, turnover.param.anc=turnover.param.anc, turnover.param.indep=turnover.param.indep, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, turnover.sigma.anc=turnover.sigma.anc, eps.param.anc=eps.param.anc, eps.param.indep=eps.param.indep, eps.sigma.indep= eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.sigma.anc=eps.sigma.anc, turn.k=turn.k, eps.k=eps.k, turnover.prob.kick=turnover.prob.kick, turnover.kick.value=turnover.kick.value, eps.prob.kick=eps.prob.kick, eps.kick.value=eps.kick.value)), NtaxGreaterThan0)))/n.sims
	}else{
		SimPhi<-function(nstart){
			sim.tree<-DeepSim(max.time=max.time, turnover.param.anc=turnover.param.anc, turnover.param.indep=turnover.param.indep, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, turnover.sigma.anc=turnover.sigma.anc, eps.param.anc=eps.param.anc, eps.param.indep=eps.param.indep, eps.sigma.indep= eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.sigma.anc=eps.sigma.anc, turn.k=turn.k, eps.k=eps.k, turnover.prob.kick=turnover.prob.kick, turnover.kick.value=turnover.kick.value, eps.prob.kick=eps.prob.kick, eps.kick.value=eps.kick.value)
			tmp<-NtaxGreaterThan0(sim.tree)
			tmp
		}
		phi.sims<-mclapply(1:n.sims, SimPhi, mc.cores=n.cores)
		actual.phi<-sum(unlist(phi.sims))/n.sims
	}
	set.seed(save.seed)
	return(actual.phi)
}

NtaxGreaterThan0<-function(sim.phylo) {
	ntax<-length(sim.phylo$tip.label)
	print(ntax)
	if(ntax>0) {
		return(TRUE)
	}
}

######################################################################################################################################
######################################################################################################################################
### This function computes Psi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Psi <- function(s, t, f, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits, t.edge) {
	#Equation broken up for ease of debugging. The extreme number of input makes it difficult to see what is going on:
	#res = first part of Morlon et al 2011 eq. 3
	res <- exp((SetDiversification(stop.time=t, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*t)-(SetDiversification(stop.time=s, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*s))
	#res2 = numerator of Morlon et al 2011 eq. 3
	res2 <- IntegrateDiversificationOverTime.int.int(stop.time=s, start.time=t, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)
	#res3 = inside brackets of Morlon et al 2011 eq. 3
	res3 <- res2/(1/f+IntegrateDiversificationOverTime.int.int(stop.time=0, start.time=s, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge))
	#res = Morlon et al 2011 eq. 3
	res <- res*(abs(1+res3)^(-2))
	return(res)
}

######################################################################################################################################
######################################################################################################################################
### Smoothing simulator function for inheritance or ratchet models  -- obtains a dataframe of quantiles for each edge:
######################################################################################################################################
######################################################################################################################################

GetQuantiles<-function(phylo){
	
	#Initiates the dataframe with the root given the 50th quantile:
	param.quantiles <- data.frame(node=Ntip(phylo)+1, turnover=0.5, eps=0.5)
	
	nbtips <- Ntip(phylo)
	for (j in 1:(nbtips-1)){
		node <- (nbtips+j)
		edges <- phylo$edge[phylo$edge[,1]==node,]
		node.anc <- phylo$edge[which(phylo$edge[,2]==node),1]
		if (length(node.anc)==0) {
			node.anc <- node
			#assume params of root = params of root parent
		}
		else {
			#Generates a column of randomly drawn quantiles for both parameters:
			param.quantiles <- rbind(param.quantiles, data.frame(node=node, turnover=runif(1,0,1), eps=runif(1,0,1)))
		}			
	}
	return(param.quantiles)
}

######################################################################################################################################
######################################################################################################################################
### Smoothing simulator function for inheritance or ratchet models -- rescales rates for an edge based on quantile:
######################################################################################################################################
######################################################################################################################################

GetNewAncParam<-function(stop.time, param.anc, sigma.indep, weight.anc.0, weight.anc.half, weight.logistic, trend.exponent, split.times, k, param.splits, param.sigma.anc, quantile.node, prob.kick, kick.value, t.edge){
	do.kick<-(prob.kick > runif(1,0,1))
	result <- (1-do.kick) * qlnorm(quantile.node,log(SetParameter(stop.time=stop.time, param.anc=param.anc, sigma.indep=sigma.indep, weight.anc.0=weight.anc.0, weight.anc.half=weight.anc.half, weight.logistic=weight.logistic, trend.exponent=trend.exponent, split.times=split.times, k=k, param.splits=param.splits, t.edge=t.edge)), param.sigma.anc) + do.kick * kick.value
	return(result)
}

######################################################################################################################################
######################################################################################################################################
### Obtains the likelihood:
######################################################################################################################################
######################################################################################################################################

GetLikelihood <- function(phylo,tot_time,f, turnover.param.root, turnover.param.indep, turnover.sigma.indep, turnover.weight.anc.0, turnover.weight.anc.half, turnover.weight.logistic, turnover.trend.exponent, turnover.sigma.anc, eps.param.root, eps.param.indep, eps.sigma.indep, eps.weight.anc.0, eps.weight.anc.half, eps.weight.logistic, eps.trend.exponent, eps.sigma.anc, split.times, turn.k, eps.k, turnover.prob.kick, turnover.kick.value, eps.prob.kick, eps.kick.value, scaled.set){

	scaled.set<-as.data.frame(scaled.set)
	indLikelihood <- c()
	nbfinal <- length(phylo)
	nbtips_tot <- 0

	turnover.param.anc <- turnover.param.root
	eps.param.anc <- eps.param.root
	ancestral.params <- data.frame(node=Ntip(phylo)+1, turnover=turnover.param.anc, eps=eps.param.anc)

	###########################REMAINING ISSUE###########################
	#Need to figure out the best way to set up timeslice. Basically, it should be an argument that is passed, as to how many and when,
	#or if any NAs, then estimate the most likely time(s) of the slice(s). Should be easyish...
	turnover.splits <- exp(rnorm(length(split.times), log(turnover.param.indep), turnover.sigma.indep))
	eps.splits <- exp(rnorm(length(split.times), log(eps.param.indep), eps.sigma.indep))
	#####################################################################
	
	from_past <- cbind(phylo$edge,node.age(phylo)$ages)
	nbtips <- Ntip(phylo)
	nbtips_tot <- nbtips_tot+nbtips
	ages <- rbind(from_past[,2:3],c(nbtips+1,0))
	ages <- ages[order(ages[,1]),]
	age <- max(ages[,2])
	#Here tj is the age of the ancestral split or the start times; sj1 is the stop time or the interval of 
	#daughter j1; and sj2 is the stop time or interval of daughter j2 until the next split. 
	#root == age.
	for (j in 1:(nbtips-1)){
		node <- (nbtips+j)
		edges <- phylo$edge[phylo$edge[,1]==node,]
		tj <- age-ages[edges[1,1],2]
		sj1 <- age-ages[edges[1,2],2]
		sj2 <- age-ages[edges[2,2],2]
		node.anc <- phylo$edge[which(phylo$edge[,2]==node),1]
		if (length(node.anc)==0) {
			node.anc <- node
			#Just needs to be non-zero so that it does not throw an error:
			t.edge=0
			#assume params of root = params of root parent
		}
		else {
			#Here stop.time is the midpoint of the last interval on ancestral branch
			last.interval <- min(split.times[which(split.times>tj)])
			stop.time <- mean(last.interval, tj)
			#Gets the length of the subtending branch to be used for any of the inheritance models:
			t.edge<-phylo$edge.length[which(phylo$edge[,2]==node)]
			#Here we want to just scale the tree based on the input list and a given inheritance parameter:
			ancestral.params <- rbind(ancestral.params, data.frame(node=node, turnover=GetNewAncParam(stop.time, param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$turnover, sigma.indep=turnover.sigma.indep, weight.anc.0=turnover.weight.anc.0, weight.anc.half=turnover.weight.anc.half,weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits, param.sigma.anc=turnover.sigma.anc, quantile.node=scaled.set[which(scaled.set$node==node),]$turnover, prob.kick=turnover.prob.kick, kick.value=turnover.kick.value,t.edge=t.edge), eps=GetNewAncParam(stop.time, param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$eps, sigma.indep=eps.sigma.indep, weight.anc.0=eps.weight.anc.0, weight.anc.half=eps.weight.anc.half, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits, param.sigma.anc=eps.sigma.anc, quantile.node=scaled.set[which(scaled.set$node==node),]$eps, prob.kick=eps.prob.kick, kick.value=eps.kick.value,t.edge=t.edge)))
		}
		indLikelihood <- c(indLikelihood,SetBirth(stop.time=tj, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*Psi(s=sj1,t=tj, f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)*Psi(s=sj2,t=tj,f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge))
	}
	#Modification of Eq. 1 from Morlon et al 2011 so the calculation is done in logspace:
	data_lik <- sum(log(indLikelihood))+(nbtips_tot*log(f))
	t.edge=0
#	Phi.value <- PhiMorlon(stop.time=tot_time,f=f, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc.0=turnover.weight.anc.0, turnover.weight.anc.half=turnover.weight.anc.half,turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc.0=eps.weight.anc.0, eps.weight.anc.half=eps.weight.anc.half, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits, t.edge=t.edge)
#	final_lik <- data_lik - log(1-Phi.value)

	return(data_lik)
}

#Print function
print.diversity<-function(x,...){
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips,x$f, row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax","SampleFrac")
	cat("\nFit\n")
	print(output)
	cat("\n")
	param.est <- matrix(0,10,2)
	param.est[,1] <- x$solution[c(1,2,7,8,9,13,17,18,15,5)]
	param.est[,2] <- x$solution[c(3,4,10,11,12,14,19,20,16,6)]
	param.est <- data.frame(param.est, row.names=c("rate","rate.indep","weight.anc.0","weight.anc.half","weight.anc.sigma","time.slice.sigma","ratchet.prob","ratchet.value","exponent","K"))
	names(param.est) <- c("turnover", "extinct.frac")
	cat("Rates\n")
	print(param.est)
	cat("\n")
}

######################################################################################################################################
######################################################################################################################################
### Code used for testing purposes:
######################################################################################################################################
######################################################################################################################################

#phy<-rcoal(130)
#phy<-read.tree("~/Dropbox/CollabBeaulieu/GeneralDiversification/test.tre")
#phy<-drop.tip(phy, paste("t", c(1:7), sep=""))
#phy<-drop.tip(phy, paste("t", c(1:5), sep=""))
#phy<-tree
#f=1
#Rprof()
#tree.set<-GetQuantiles(phylo=phy)
#scaled.set<-tree.set
#scaled.set[,2]<-qlnorm(tree.set[,2], log(.85), 0.1)
#scaled.set[,3]<-qlnorm(tree.set[,3], log(.25), 0.00)
#us.DD.LH <- (GetLikelihood(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=0.5, turnover.sigma.indep=0,turnover.sigma.anc=0.05, turnover.weight.anc.0=.75, turnover.weight.logistic=0, turnover.trend.exponent=0, eps.param.indep=0.70, eps.sigma.indep=0, eps.weight.anc.0=0, eps.weight.logistic=0, eps.trend.exponent=0, eps.sigma.anc=0, split.times=branching.times(phy), turn.k=1, eps.k=1, scaled.set=scaled.set))
#Rprof(NULL)
#summaryRprof()
#0.4948420525 0.7078933716
						  
#phy<-tree
#turnover.sigma.anc<-0:10*0.1
#turnover.weight.anc.0<-0:10*0.1
#quantile.set<-lapply(1:100,GetQuantileSet)
#tmp<-numeric(100)
#loglik<-c()
#weight.anc.0<-c()
#sigma.anc<-c()
#for(k in 1:2){
#	for(i in 1:11){
#		for(j in 1:100){
#			tmp[j]<-GetLikelihood(phylo=phy,tot_time=max(branching.times(phy)),f=88/120, turnover.param.indep=0.539856, turnover.sigma.indep=0,turnover.sigma.anc=turnover.sigma.anc[k], turnover.weight.anc.0=turnover.weight.anc.0[i], turnover.weight.logistic=0, turnover.trend.exponent=0, eps.param.indep=0.8574219, eps.sigma.indep=0, eps.weight.anc.0=0, eps.weight.logistic=0, eps.trend.exponent=0, eps.sigma.anc=0, split.times=branching.times(phy), turn.k=1, eps.k=1, scaled.set=quantile.set[j])
#		}
#		loglik <-c(loglik,mean(tmp))
#		weight.anc.0<-c(weight.anc.0,turnover.weight.anc.0[i])
#		sigma.anc<-c(sigma.anc,turnover.sigma.anc[k])
#	}
#}




