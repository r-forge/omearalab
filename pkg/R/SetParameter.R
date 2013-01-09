#Overall goal is to have turnover rate or extinction fraction vary as a result of various
#   parameters.
#Note that with logistic growth, it makes the most sense to pass in split.times as actual
#   split.times (including those with no descendants, but in practice, split.times reconstructed
#   on the phylogeny will probably be used instead. Rabosky and Lovette 2008, in a response
#   argue this is fairly okay (the two measures are correlated).
#Inputs:
#phy =the phylogeny
#f = sampling fraction
#model = whether a yule type model is to be used or birth-death
#logistic = whether or not the logistic growth model is to be specified
#inherit = whether rates are to be inherited
#brown = whether rates vary according to a BM process
#exp = whether rates vary exponentially

library(picante)
library(ape)
library(nloptr)
library(multicore)

GeneralDiversity<-function(phy, f=1, model=c("yule", "bd"), turnover.logistic=FALSE, eps.logistic=FALSE, turnover.exp=FALSE, eps.exp=FALSE, turnover.inherit=FALSE, eps.inherit=FALSE, turnover.slice=FALSE, eps.slice=FALSE, n.cores=NULL) {
	
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
		}
		if(model=="bd"){
			turnover.param = TRUE
			eps.param = TRUE
		}
	}
	#Sets two important inputs:
	split.times=sort(branching.times(phy), decreasing=TRUE)
	tot_time=max(branching.times(phy))
	#Makes a vector of parameters that are going to be estimated:
	pars=c(turnover.param, turnover.inherit, eps.param, eps.inherit, turnover.logistic, eps.logistic, turnover.inherit, turnover.inherit, eps.inherit, eps.inherit, turnover.slice, eps.slice, turnover.exp, eps.exp)
	#An index for use later:
	tmp<-pars
	#The trues specify the number of parameters:
	pars[pars==T] <- 1:length(pars[pars==T])
	np<-max(pars)
	#All parameters not estimated are set to 1+max number of estimated parameters:
	pars[pars==0]<-max(pars)+1
	#Function used for optimizing parameters:
	DevOptimize <- function(p, pars, phy, tot_time, f, turnover.inherit, turnover.weight.logistic, eps.inherit, eps.weight.logistic, split.times, quantile.set, n.cores) {
		#Generates the final vector with the appropriate parameter estimates in the right place:
		model.vec <- numeric(length(pars))
		model.vec[] <- c(p, 0)[pars]
		#Now set entries in model.vec to the defaults if we are not estimating them:
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
		if(pars[8] < max(pars) | pars[10] < max(pars) | pars[11] < max(pars) | pars[12] < max(pars)){
			if(is.null(n.cores)){
				tot.logl <-c()
				for(i in 1:length(quantile.set)){
					tot.logl <- c(tot.logl,GetLikelihood(phylo=phy, tot_time, f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[11], turnover.weight.anc=model.vec[7], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[13], turnover.sigma.anc=model.vec[8], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[12], eps.weight.anc=model.vec[9], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[14], eps.sigma.anc=model.vec[10], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], scaled.set=quantile.set[[i]]))
				}
				logl<-mean(tot.logl[is.finite(tot.logl)],na.rm=T)
				if(!is.finite(logl)){
					return(10000000)
				}
			}
			else{
				SimSigma<-function(nstart){
					tmp = matrix(,1,ncol=1)
					tot.logl <- GetLikelihood(phylo=phy, tot_time, f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[11], turnover.weight.anc=model.vec[7], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[13], turnover.sigma.anc=model.vec[8], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[12], eps.weight.anc=model.vec[9], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[14], eps.sigma.anc=model.vec[10], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], scaled.set=quantile.set[[nstart]])
					tmp[,1] = tot.logl
					tmp
				}
				sim.set<-mclapply(1:length(quantile.set), SimSigma, mc.cores=n.cores)
				tot.logl<-unlist(sim.set)
				logl<-mean(tot.logl[is.finite(tot.logl)],na.rm=T)
			}
			return(-logl)
		}
		else{
			logl <- GetLikelihood(phylo=phy, tot_time, f, turnover.param.root=model.vec[1], turnover.param.indep=model.vec[2], turnover.sigma.indep=model.vec[11], turnover.weight.anc=model.vec[7], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=model.vec[13], turnover.sigma.anc=model.vec[8], eps.param.root=model.vec[3], eps.param.indep=model.vec[4], eps.sigma.indep=model.vec[12], eps.weight.anc=model.vec[9], eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=model.vec[14], eps.sigma.anc=model.vec[10], split.times=split.times, turn.k=model.vec[5], eps.k=model.vec[6], scaled.set=quantile.set)
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
	init.pars=c(turnover.param, FALSE, eps.param, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
	init.tmp<-init.pars
	init.pars[init.pars==T] <- 1:length(init.pars[init.pars==T])
	init.pars[init.pars==0]<-max(init.pars)+1
	init.set.pars <- c(1,1,0.5,0.5,Ntip(phy)*2,Ntip(phy)*2,0.5,0.5,0.5,0.5,0.5,0.5,0,0)
	init.ip <- init.set.pars[init.tmp==TRUE]
	init.set.lower <- c(0,0,0,0,-1e6,-1e6,0,0,0,0,0,0,-10,-10)
	lower.init <- init.set.lower[init.tmp==TRUE]
	init.set.upper <- c(100,100,100,100,1e6,1e6,1,10,1,10,10,10,10,10)
	upper.init <- init.set.upper[init.tmp==TRUE]
	init = nloptr(x0=init.ip, eval_f=DevOptimize, lb=lower.init, ub=upper.init, opts=opts, pars=init.pars, phy=phy, tot_time=tot_time, f=f, turnover.inherit=FALSE, turnover.weight.logistic=FALSE, eps.inherit=FALSE, eps.weight.logistic=FALSE, split.times=split.times, quantile.set=init.quantile.set, n.cores=n.cores)

	cat("Finished. Begin thorough search...", "\n")
	
	#Now begin more thorough search using estimates from constant bd model:
	def.set.pars <- c(init$solution[1],init$solution[1],init$solution[2],init$solution[2],Ntip(phy)*2,Ntip(phy)*2,0.5,0.5,0.5,0.5,0.5,0.5,0,0)
	ip <- def.set.pars[tmp==TRUE]
	def.set.lower <- c(0,0,0,0,-1e6,-1e6,0,0,0,0,0,0,-10,-10)
	lower <- def.set.lower[tmp==TRUE]
	def.set.upper <- c(100,100,100,100,1e6,1e6,1,10,1,10,10,10,10,10)
	upper <- def.set.upper[tmp==TRUE]

	#If inheritance model is chosen, create a list of quantiles for each edge in the tree:
	if(turnover.inherit==TRUE | eps.inherit==TRUE){ 
		GetQuantileSet<-function(nstarts){
			GetQuantiles(phylo=phy)
		}
		quantile.set<-lapply(1:1000, GetQuantileSet)
	}
	else{
		quantile.set<-GetQuantiles(phylo=phy)
	}
	
	out = nloptr(x0=ip, eval_f=DevOptimize, lb=lower, ub=upper, opts=opts, pars=pars, phy=phy, tot_time=tot_time, f=f, turnover.inherit=turnover.inherit, turnover.weight.logistic=turnover.logistic, eps.inherit=eps.inherit, eps.weight.logistic=eps.logistic, split.times=split.times, quantile.set=quantile.set, n.cores=n.cores)
	loglik <- -out$objective
	#Recreate the model vector for use in the print function:
	solution <- numeric(length(pars))
	solution[] <- c(out$solution, 0)[pars]

	obj = list(loglik = loglik, AIC = -2*loglik+2*np,AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par = pars, opts=opts, f=f, phy=phy, iterations=out$iterations) 
	
	class(obj)<-"diversity"		
	return(obj)	
}

SetParameter <- function(stop.time, param.anc, sigma.indep, weight.anc, weight.logistic, trend.exponent, split.times, k, param.splits) {
	position <- max(which(split.times>=stop.time),1)
	n.taxa <- 1+length(which(split.times>=stop.time))
	param.indep <- param.splits[position]
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc) #same as autoregressive model
	logistic.scaling <- 1 - (weight.logistic * (n.taxa / k))
	param.mean <- (param.starting * (n.taxa ^ trend.exponent)) * logistic.scaling
	return(param.mean)
}

SetBirth <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	turnover <- SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits)
	extinction.fraction <- SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits)
	return(SetBirthGivenParameters(turnover, extinction.fraction))
}

SetBirthGivenParameters <- function(turnover, extinction.fraction) {
	birth.expected <- turnover / (1 + extinction.fraction)
	return(birth.expected)
}

SetDeath <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	turnover<-SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits)
	extinction.fraction<-SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits)
	return(SetDeathGivenParameters(turnover, extinction.fraction))
}

SetDeathGivenParameters <- function(turnover, extinction.fraction) {
	death.expected<-(extinction.fraction * turnover) / (1 + extinction.fraction)
	return(death.expected)
}

SetDiversification <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	turnover<-SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits)
	extinction.fraction<-SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits)
	return(SetBirthGivenParameters(turnover, extinction.fraction) - SetDeathGivenParameters(turnover, extinction.fraction))
}

######################################################################################################################################
######################################################################################################################################
### This is our solution for integating diversification  
######################################################################################################################################
######################################################################################################################################

IntegrateSetDiversificationExpected <- function(lower, upper, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic,turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {	
	if (turnover.weight.logistic>0 | eps.weight.logistic>0 | turnover.trend.exponent>0 | eps.trend.exponent>0) {
		boundaries<-split.times[which(split.times>lower)]
		boundaries<-boundaries[which(boundaries<upper)]
		lower.bounds<-sort(c(lower, boundaries))
		upper.bounds<-sort(c(boundaries, upper))
		integration.result<-0
		for (i in sequence(length(lower.bounds))) {
			#piecewise equation: (exp(b-d)*t-exp(b-d)*s)*b / (b-d) 
			integration.result <- integration.result + ((exp(SetDiversification(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*upper.bounds[i])-exp(SetDiversification(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*lower.bounds[i]))*SetBirth(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits))/SetDiversification(stop.time=(lower.bounds[i]+upper.bounds[i])/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
		}
		return(integration.result)
	}
	else {
		#piecewise equation: (exp(b-d)*t-exp(b-d)*s)*b / (b-d)
		integration.result<-((exp(SetDiversification(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*upper)-exp(SetDiversification(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*lower))*SetBirth(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)) / SetDiversification(stop.time=(lower+upper)/2, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	}
	return(integration.result)
}

IntegrateDiversificationOverTime <- function(stop.time, start.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	return(IntegrateSetDiversificationExpected(lower=stop.time, upper=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits))
}

IntegrateDiversificationOverTime.int.int <- function(stop.time, start.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	result<-IntegrateDiversificationOverTime(stop.time=stop.time, start.time=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	return(result)
}

######################################################################################################################################
######################################################################################################################################
### This function computes Phi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Phi <- function(stop.time, f, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	res <- (1-exp((SetDiversification(stop.time=stop.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*stop.time)-(SetDiversification(stop.time=0, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*0)) / (1/f+IntegrateDiversificationOverTime.int.int(stop.time=0,start.time=stop.time, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)))
	return(res)	
}

######################################################################################################################################
######################################################################################################################################
### This function computes Psi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Psi <- function(s, t, f, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	#Equation broken up for ease of debugging. The extreme number of input makes it difficult to see what is going on:
	#res = first part of Morlon et al 2011 eq. 3
	res <- exp((SetDiversification(stop.time=t, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*t)-(SetDiversification(stop.time=s, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*s))
	#res2 = numerator of Morlon et al 2011 eq. 3
	res2 <- IntegrateDiversificationOverTime.int.int(stop.time=s, start.time=t, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	#res3 = inside brackets of Morlon et al 2011 eq. 3
	res3 <- res2/(1/f+IntegrateDiversificationOverTime.int.int(stop.time=0, start.time=s, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	#res = Morlon et al 2011 eq. 3
	res <- res*(abs(1+res3)^(-2))
	return(res)
}

######################################################################################################################################
######################################################################################################################################
### For inheritance model obtains a dataframe of quantiles for each edge:
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
### For inheritance model rescales tree based on weight of the ancestor parameter:
######################################################################################################################################
######################################################################################################################################

GetNewAncParam<-function(stop.time, param.anc, sigma.indep, weight.anc, weight.logistic, trend.exponent, split.times, k, param.splits, param.sigma.anc, quantile.node){
	result<-qlnorm(quantile.node,log(SetParameter(stop.time=stop.time, param.anc=param.anc, sigma.indep=sigma.indep, weight.anc=weight.anc, weight.logistic=weight.logistic, trend.exponent=trend.exponent, split.times=split.times, k=k, param.splits=param.splits)), param.sigma.anc)
	return(result)
}

#GetNewAncParam<-function(stop.time, param.anc, sigma.indep, weight.anc, weight.logistic, trend.exponent, split.times, k, param.splits, param.sigma.anc){
#	result<-exp(rnorm(1,log(SetParameter(stop.time=stop.time, param.anc=param.anc, sigma.indep=sigma.indep, weight.anc=weight.anc, weight.logistic=weight.logistic, trend.exponent=trend.exponent, split.times=split.times, k=k, param.splits=param.splits)), param.sigma.anc))
#	return(result)
#}

######################################################################################################################################
######################################################################################################################################
### Obtains the likelihood:
######################################################################################################################################
######################################################################################################################################

GetLikelihood <- function(phylo,tot_time,f, turnover.param.root, turnover.param.indep, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, turnover.sigma.anc, eps.param.root, eps.param.indep, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, eps.sigma.anc, split.times, turn.k, eps.k, scaled.set){

	scaled.set<-as.data.frame(scaled.set)
	indLikelihood <- c()
	nbfinal <- length(phylo)
	nbtips_tot <- 0

	turnover.param.anc <- turnover.param.root
	eps.param.anc <- eps.param.root
	ancestral.params <- data.frame(node=Ntip(phylo)+1, turnover=turnover.param.anc,eps=eps.param.anc)

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
			#assume params of root = params of root parent
		}
		else {
			#Here stop.time is the midpoint of the last interval on ancestral branch
			last.interval <- min(split.times[which(split.times>tj)])
			stop.time <- mean(last.interval, tj)
			#Here we want to just scale the tree based on the input list and a given inheritance parameter:
			ancestral.params <- rbind(ancestral.params, data.frame(node=node, turnover=GetNewAncParam(stop.time, param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$turnover, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits, param.sigma.anc=turnover.sigma.anc, quantile.node=scaled.set[which(scaled.set$node==node),]$turnover), eps=GetNewAncParam(stop.time, param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$eps, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits, param.sigma.anc=eps.sigma.anc, quantile.node=scaled.set[which(scaled.set$node==node),]$eps)))
		}
		indLikelihood <- c(indLikelihood,SetBirth(stop.time=tj, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*Psi(s=sj1,t=tj, f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)*Psi(s=sj2,t=tj,f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	}
	#Modification of Eq. 1 from Morlon et al 2011 so the calculation is done in logspace:
	data_lik <- sum(log(indLikelihood))+(nbtips_tot*log(f))
	Phi <- Phi(stop.time=tot_time,f=f, turnover.param.anc=turnover.param.anc, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	final_lik <- data_lik - log(1-Phi)

	return(final_lik)
}

#Print function
print.diversity<-function(x,...){
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips,x$f, row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax","SampleFrac.")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	param.est<- matrix(0,7,2)
	param.est[,1] <- x$solution[c(1,2,7,8,11,13,5)]
	param.est[,2] <- x$solution[c(3,4,8,9,12,14,6)]
	param.est <- data.frame(param.est, row.names=c("rate","rate.indep","weight.anc","weight.anc.sigma","time.slice.sigma","exponent","K"))
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
#us.DD.LH <- (GetLikelihood(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=0.5, turnover.sigma.indep=0,turnover.sigma.anc=0.05, turnover.weight.anc=.75, turnover.weight.logistic=0, turnover.trend.exponent=0, eps.param.indep=0.70, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.exponent=0, eps.sigma.anc=0, split.times=branching.times(phy), turn.k=1, eps.k=1, scaled.set=scaled.set))
#Rprof(NULL)
#summaryRprof()
#0.4948420525 0.7078933716
						  
#phy<-tree
#turnover.sigma.anc<-0:10*0.1
#turnover.weight.anc<-0:10*0.1
#quantile.set<-lapply(1:100,GetQuantileSet)
#tmp<-numeric(100)
#loglik<-c()
#weight.anc<-c()
#sigma.anc<-c()
#for(k in 1:2){
#	for(i in 1:11){
#		for(j in 1:100){
#			tmp[j]<-GetLikelihood(phylo=phy,tot_time=max(branching.times(phy)),f=88/120, turnover.param.indep=0.539856, turnover.sigma.indep=0,turnover.sigma.anc=turnover.sigma.anc[k], turnover.weight.anc=turnover.weight.anc[i], turnover.weight.logistic=0, turnover.trend.exponent=0, eps.param.indep=0.8574219, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.exponent=0, eps.sigma.anc=0, split.times=branching.times(phy), turn.k=1, eps.k=1, scaled.set=quantile.set[j])
#		}
#		loglik <-c(loglik,mean(tmp))
#		weight.anc<-c(weight.anc,turnover.weight.anc[i])
#		sigma.anc<-c(sigma.anc,turnover.sigma.anc[k])
#	}
#}




