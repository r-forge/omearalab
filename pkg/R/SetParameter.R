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
#brown.trend = whether rates vary according to a BM process with a trend

####### TO DO #######

#1. Remove all the sigma.time variables from all functions
#2. Test
#3. Add in identifiability tests?
#4. Need better name

#####################

library(picante)
library(ape)
library(nloptr)

GeneralDiversity<-function(phy, f=1, model=c("yule", "bd"), turnover.logistic=TRUE, eps.logistic=TRUE, turnover.inherit=TRUE, eps.inherit=TRUE, turnover.brown=TRUE, eps.brown=TRUE, turnover.brown.trend=TRUE, eps.brown.trend=TRUE) {
	
	obj <- NULL
	#Sets the main parameters to be used in the model:
	if (is.character(model)) {
		if(model=="yule"){
			turnover.indep = TRUE
			eps.param.indep = FALSE
			eps.logistic=FALSE
			eps.inherit=FALSE
			eps.brown=FALSE
			eps.brown.trend=FALSE
		}
		if(model=="bd"){
			turnover.indep = TRUE
			eps.param.indep = TRUE
		}
	}
	#Sets two important inputs:
	split.times=branching.times(phy)
	tot_time=max(branching.times(phy))
	#Are we going to be using a logistic model?
	if(turnover.weight.logistic == TRUE | eps.weight.logistic == TRUE){
		logistic.model <- TRUE
	}
	else{
		logistic.model <- FALSE
	}	
	#Makes a vector of parameters that are going to be estimated:
	pars=c(turnover.indep, eps.param.indep, logistic.model, turnover.weight.anc, eps.weight.anc, turnover.sigma.indep, eps.sigma.indep)
	#Tack on the argument for whether the BM has a trend -- in testing it did not like the use of else so a ton of ifs here:
	if(turnover.trend==TRUE) {
		pars<-c(pars,TRUE,TRUE)
	}
	if(turnover.trend==FALSE) {
		pars<-c(pars,FALSE,FALSE)
	}
	if(eps.trend==TRUE) {
		pars<-c(pars,TRUE,TRUE)
	}
	if(eps.trend==FALSE) {
		pars<-c(pars,FALSE,FALSE)
	}
	#The trues specify the number of parameters:
	pars[pars==T]<-1:length(pars[pars==T])
	np<-max(pars)
	#All parameters not estimated are set to 1+max number of estimated parameters:
	pars[pars==0]<-max(pars)+1
	#Function used for 
	DevOptimize <- function(p, pars, phylo, tot_time, f, turnover.weight.logistic, eps.weight.logistic, split.times) {
		#Generates the final vector with the appropriate parameter estimates in the right place:
		model.vec <- numeric(length(pars))
		model.vec[] <- c(p, 0)[pars]
		#Now set entries in model.vec to the defaults if we are not estimating them:
		if(model.vec[3]==0) {
			model.vec[3] = 1
		}		
		if(model.vec[9]==0) {
			model.vec[9] = 1
		}
		if(model.vec[11]==0) {
			model.vec[11] = 1
		}
		logl <- getLikelihood.gen.bd(phylo=phy, tot_time, f, turnover.param.indep=model.vec[1], turnover.sigma.indep=model.vec[6], turnover.weight.anc=model.vec[4], turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=model.vec[8], turnover.trend.exponent=model.vec[9], eps.param.indep=model.vec[2], eps.sigma.indep=model.vec[7], eps.weight.anc=model.vec[5], eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=model.vec[10], eps.trend.exponent=model.vec[11], split.times=split.times, k=model.vec[3])
		return(-logl)
	}
	
	#ip
	#lb
	#ub
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	out = nloptr(x0=rep(ip, length.out = np), eval_f=DevOptimize, opts=opts, pars=pars, phylo=phy, tot_time=tot_time, f=f, turnover.weight.logistic=turnover.weight.logistic, eps.weight.logistic=eps.weight.logistic, split.times=split.times)
	loglik <- -out$objective
	#Recreate the model vector for use in the print function:
	solution <- numeric(length(pars))
	solution[] <- c(p, 0)[pars]
	if(model.vec[3]==0) {
		model.vec[3] = 1
	}		
	if(model.vec[9]==0) {
		model.vec[9] = 1
	}
	if(model.vec[11]==0) {
		model.vec[11] = 1
	}
	obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=2*loglik+2*np+(2*np*(np+1))/(1-np-1), solution=solution, opts=opts, f=f, phy=phy, lb=lower, ub=upper, iterations=out$iterations) 

	class(obj)<-"diversity"		
	return(obj)	
}	

SetParameter <- function(stop.time, param.anc, sigma.time, sigma.indep, weight.anc, weight.logistic, trend.scaling=0, trend.exponent=1, split.times, k, param.splits) {
	position <- max(which(split.times>=stop.time),1)
	n.taxa<-1+position
	param.indep <- param.splits[position]
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc)
	param.mean <- (param.starting + trend.scaling * (stop.time ^ trend.exponent)) * (1 - (weight.logistic * n.taxa / k))
	return(param.mean)
} 

SetBirth <- function(stop.time, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	turnover <- SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, param.splits=turnover.splits)
	extinction.fraction <- SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k, param.splits=eps.splits)
	return(SetBirthGivenParameters(turnover, extinction.fraction))
}

SetBirthGivenParameters <- function(turnover, extinction.fraction) {
	birth.expected <- turnover / (1 + extinction.fraction)
	return(birth.expected)
}

IndefiniteIntegralForDiversification <- function(time, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	
	position<-max(which(split.times>=time),1)
	n.taxa<-1+position
	turnover.param.indep <- turnover.splits[position]
	eps.param.indep <- eps.splits[position]
	
	turnover.log <- 1 - (turnover.weight.logistic * n.taxa / k)
	turnover.indep <- turnover.param.indep * (1 - turnover.weight.anc)
	turnover.anc <- turnover.param.anc * turnover.weight.anc
	eps.log <- 1 - (eps.weight.logistic * n.taxa / k)
	eps.indep <- eps.param.indep * (1 - eps.weight.anc)
	eps.anc <- eps.param.anc * eps.weight.anc
	
	turnover <- turnover.log*(((turnover.trend.scaling*(time^(turnover.trend.exponent+1)))/(turnover.trend.exponent+1))+(turnover.anc*time)+(turnover.indep*time))
	extinction.fraction <- eps.log*(((eps.trend.scaling*(time^(eps.trend.exponent+1)))/(eps.trend.exponent+1))+(eps.anc*time)+(eps.indep*time))
	diversification <- (turnover / (1 + extinction.fraction)) - ((extinction.fraction * turnover) / (1 + extinction.fraction))	
	
	return(diversification)
}

IntegrateDiversificationOverTime.int.0 <- function(start.time, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	return(exp(IndefiniteIntegralForDiversification(time=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))*SetBirth(stop.time=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
}

IntegrateDiversificationOverTime.int.int <- function(stop.time, start.time, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	return(integrate(Vectorize(IntegrateDiversificationOverTime.int.0, "start.time"), lower=stop.time, upper=start.time, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits, stop.on.error=FALSE)$value)
}

######################################################################################################################################
######################################################################################################################################
### This function computes Phi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Phi <- function(stop.time, f, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	res <- (1-exp(IndefiniteIntegralForDiversification(time=(stop.time-0), turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)) / (1/f+IntegrateDiversificationOverTime.int.int(stop.time=0,start.time=stop.time, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)))
	return(res)	
}

######################################################################################################################################
######################################################################################################################################
### This function computes Psi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Psi <- function(s, t, f, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	#Equation broken up for ease of debugging. The extreme number of input makes it difficult to see what is going on:
	#res = first part of Morlon et al 2011 eq. 3
	res <- exp(IndefiniteIntegralForDiversification(time=(t-s), turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	#res2 = numerator of Morlon et al 2011 eq. 3
	res2 <- IntegrateDiversificationOverTime.int.int(stop.time=s, start.time=t, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	#res3 = inside brackets of Morlon et al 2011 eq. 3
	res3 <- res2/(1/f+IntegrateDiversificationOverTime.int.int(stop.time=0, start.time=s, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	#res = Morlon et al 2011 eq. 3
	res <- res*(abs(1+res3)^(-2))
	return(res)
}


getGeometricMeanParam <- function(start.time, stop.time, param.anc, sigma.time, sigma.indep, weight.anc, weight.logistic, trend.scaling=0, trend.exponent=1, split.times, k, param.splits) {
	ancestral.rate.geometric.mean <- exp(integrate(Vectorize(SetParameter, "stop.time"), lower=stop.time, upper=start.time, param.anc=param.anc, sigma.time=sigma.time, sigma.indep=sigma.indep, weight.anc=weight.anc, weight.logistic=weight.logistic, trend.scaling=trend.scaling, trend.exponent=trend.exponent, split.times=split.times, k=k, param.splits=param.splits)$value)
	return(ancestral.rate.geometric.mean)
}

######################################################################################################################################
######################################################################################################################################
### Obtains the likelihood:
######################################################################################################################################
######################################################################################################################################

getLikelihood.gen.bd <- function(phylo,tot_time,f, turnover.param.indep, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k){
	indLikelihood <- c()
	nbfinal <- length(phylo)
	nbtips_tot <- 0
	
	turnover.param.anc <- turnover.param.indep
	eps.param.anc <- eps.param.indep
	#at each step, take the geometric mean of the rate on the edge. This is the "state" of the edge. The param.anc at the edges descendant node is this mean with edge.length * sigma.time 
	#the .sigma.indep only get called here. Could clean up the inputs by removing them from all other steps.
	ancestral.params <- data.frame(node=Ntip(phylo)+1, turnover=turnover.param.anc, eps=eps.param.anc)
	turnover.splits <- abs(rnorm(length(split.times), 0, turnover.sigma.indep)+turnover.param.indep)
	eps.splits <- abs(rnorm(length(split.times), 0, eps.sigma.indep)+eps.param.indep)
	#
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
		node.anc <- phy$edge[which(phy$edge[,2]==node),1]
		if (length(node.anc)==0) {
			node.anc <- node
		#assume params of root = params of root parent
		}
		else {
			ancestral.params <- rbind(ancestral.params, data.frame(node=node, 
																 turnover=getGeometricMeanParam(start.time=tj, stop.time=tj+phy$edge.length[j], param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$turnover, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k, param.splits=turnover.splits), 
																 eps=getGeometricMeanParam(start.time=tj, stop.time=tj+phy$edge.length[j], param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$eps, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k, param.splits=eps.splits)))
		}
		indLikelihood <- c(indLikelihood,SetBirth(stop.time=tj, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)*Psi(s=sj1,t=tj, f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)*Psi(s=sj2,t=tj,f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	}
	#Can be removed -- we do not care about stem age
	#indLikelihood<-c(indLikelihood,Psi(s=age,t=tot_time,f=f,turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k))
	#Eq. 1 from Morlon et al 2011:
	data_lik <- prod(indLikelihood)*f^nbtips_tot
	Phi <- Phi(stop.time=tot_time,f=f, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	final_lik <- data_lik/(1-Phi)
	
	return(log(final_lik))
}

#Print function
print.diversity<-function(x,...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips,f, row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax","SampleFrac.")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	cat("Carrying capacity")
	if(x$solution==1) {
		print("logistic model turned off")
	}
	else{
		print(x$solution)
	}
	cat("\n")
	param.est<- matrix(,5,2)
	param.est[,1] <- x$solution[c(1,4,6,8,9)]
	param.est[,2] <- x$solution[c(2,5,7,10,11)]
	param.est <- data.frame(param.est, row.names="rate", "weight.anc", "sigma", "trend.scaling", "trend.exponent")
	names(param.est) <- c("turnover", "extinction.fraction")
	cat("Rates\n")
	print(param.est)
	cat("\n")
}




######################################################################################################################################
######################################################################################################################################
### Code used for testing purposes:
######################################################################################################################################
######################################################################################################################################

#phy<-rcoal(1000)
phy<-read.tree("~/Dropbox/CollabBeaulieu/GeneralDiversification/test.tre")
#phy<-drop.tip(phy, paste("t", c(1:7), sep=""))
#phy<-drop.tip(phy, paste("t", c(1:5), sep=""))

us.DD.LH <- (getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=0.5, turnover.sigma.time=0.0, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=0, eps.sigma.time=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, split.times=branching.times(phy), k=100))
us.bd.LH <- (getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=0.5, turnover.sigma.time=0, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=0, eps.sigma.time=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, split.times=branching.times(phy), k=1))

print(paste("us.bd.LH", us.bd.LH))
print(paste("us.DD.LH", us.DD.LH))
print(paste("us.DD.LH - us.bd.LH", us.DD.LH - us.bd.LH))
