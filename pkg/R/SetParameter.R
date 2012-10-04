#Overall goal is to have turnover rate or extinction fraction vary as a result of various
#   parameters.
#Note that with logistic growth, it makes the most sense to pass in split.times as actual
#   split.times (including those with no descendants, but in practice, split.times reconstructed
#   on the phylogeny will probably be used instead. Rabosky and Lovette 2008, in a response
#   argue this is fairly okay (the two measures are correlated).
#Inputs:
#	t = duration of time
#	param.root = starting parameter value at the root
#	param.anc = starting parameter value for start of interval
#	sigma.param = variance parameter
#	weight.anc = determines rate inheritance
#	weight.logistic = decides whether the model is logistic
#	trend.scaling = 
#	trend.exponent = 
#	split.times = number of taxa at a start time
#	k = carrying capacity

#general.diversity<-function(, params.turnover.constant, params.eps.constant){

#Inputs to getLikelihood...:(phylo,tot_time,f, turnover.param.indep, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k)

#init<-c(lamb_par,mu_par,f)
#p<-length(init)    
#dev.lh <- function(init) {
#	lamb_par <- init[1:length(lamb_par)]
#	mu_par <- init[(1+length(lamb_par)):(length(init)-1)]
#	f<-init[length(init)]
#	f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
#	f.mu.par<-function(x){abs(f.mu(x,mu_par))}
#	LH <- getLikelihood.gen.bd(phylos,tot_time,f.lamb.par,f.mu.par,min(abs(f),1),cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
#	return(-LH)
#}
#
#opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
#out = nloptr(x0=rep(ip, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)
#res <-list(LH = -out$objective,aicc=2*temp$objective+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$solution[1:length(lamb_par)], mu_par=temp$solution[(1+length(lamb_par)):(length(init)-1)],f=min(abs(temp$par[length(init)]),1))
#}

#return(res)

#}	

SetParameter <- function(stop.time, param.anc, sigma.time, sigma.indep, weight.anc, weight.logistic, trend.scaling=0, trend.exponent=1, split.times, k, param.splits) {
	position<-max(which(split.times>=stop.time),1)
	n.taxa<-1+position
	param.indep <- param.splits[position]
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc)
	param.mean <- (param.starting + trend.scaling * (stop.time ^ trend.exponent)) * (1 - (weight.logistic * n.taxa / k))
	return(param.mean)
} 

SetBirth <- function(stop.time, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	turnover<-SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, param.splits=turnover.splits)
	extinction.fraction<-SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k, param.splits=eps.splits)
	return(SetBirthGivenParameters(turnover, extinction.fraction))
}

SetBirthGivenParameters <- function(turnover, extinction.fraction) {
	birth.expected<-turnover / (1 + extinction.fraction)
	return(birth.expected)
}

IndefiniteIntegralForDiversification <- function(time, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	
	position<-max(which(split.times>=time),1)
	n.taxa<-1+position
	turnover.param.indep <- turnover.splits[position]
	eps.param.indep <- eps.splits[position]
	
	turnover.log <- (1 - (turnover.weight.logistic * n.taxa / k))
	turnover.indep <- (turnover.param.indep * (1 - turnover.weight.anc))
	turnover.anc <- (turnover.param.anc * turnover.weight.anc)
	eps.log <- (1 - (eps.weight.logistic * n.taxa / k))
	eps.indep <- (eps.param.indep * (1 - eps.weight.anc))
	eps.anc <- (eps.param.anc * eps.weight.anc)
	
	turnover<-turnover.log*(((turnover.trend.scaling*(time^(turnover.trend.exponent+1)))/(turnover.trend.exponent+1))+(turnover.anc*time)+(turnover.indep*time))
	extinction.fraction<-eps.log*(((eps.trend.scaling*(time^(eps.trend.exponent+1)))/(eps.trend.exponent+1))+(eps.anc*time)+(eps.indep*time))
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

Phi<-function(stop.time, f, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	res<-(1-exp(IndefiniteIntegralForDiversification(time=(stop.time-0), turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)) / (1/f+IntegrateDiversificationOverTime.int.int(stop.time=0,start.time=stop.time, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)))
	return(res)	
}

######################################################################################################################################
######################################################################################################################################
### This function computes Psi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Psi<-function(s, t, f, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k, turnover.splits, eps.splits) {
	#Equation broken up for ease of debugging. The extreme number of input makes it difficult to see what is going on:
	#res = first part of Morlon et al 2011 eq. 3
	res<-exp(IndefiniteIntegralForDiversification(time=(t-s), turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	#res2 = numerator of Morlon et al 2011 eq. 3
	res2<-IntegrateDiversificationOverTime.int.int(stop.time=s, start.time=t, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	#res3 = inside brackets of Morlon et al 2011 eq. 3
	res3<-res2/(1/f+IntegrateDiversificationOverTime.int.int(stop.time=0, start.time=s, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	#res = Morlon et al 2011 eq. 3
	res<-res*(abs(1+res3)^(-2))
	return(res)
}


getGeometricMeanParam<-function(start.time, stop.time, param.anc, sigma.time, sigma.indep, weight.anc, weight.logistic, trend.scaling=0, trend.exponent=1, split.times, k, param.splits) {
	ancestral.rate.geometric.mean<-exp(integrate(Vectorize(SetParameter, "stop.time"), lower=stop.time, upper=start.time, param.anc=param.anc, sigma.time=sigma.time, sigma.indep=sigma.indep, weight.anc=weight.anc, weight.logistic=weight.logistic, trend.scaling=trend.scaling, trend.exponent=trend.exponent, split.times=split.times, k=k, param.splits=param.splits)$value)
	return(ancestral.rate.geometric.mean)
}

######################################################################################################################################
######################################################################################################################################
### Obtains the likelihood:
######################################################################################################################################
######################################################################################################################################

getLikelihood.gen.bd<-function(phylo,tot_time,f, turnover.param.indep, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k){
	indLikelihood<-c()
	nbfinal<-length(phylo)
	nbtips_tot<-0
	
	turnover.param.anc<-turnover.param.indep
	eps.param.anc<-eps.param.indep
	#at each step, take the geometric mean of the rate on the edge. This is the "state" of the edge. The param.anc at the edges descendant node is this mean with edge.length * sigma.time 
	#the .sigma.indep only get called here. Could clean up the inputs by removing them from all other steps.
	ancestral.params<-data.frame(node=Ntip(phylo)+1, turnover=turnover.param.anc, eps=eps.param.anc)
	turnover.splits<-abs(rnorm(length(split.times), 0,  turnover.sigma.indep)+turnover.param.indep)
	eps.splits<-abs(rnorm(length(split.times), 0,  eps.sigma.indep)+eps.param.indep)
	#
	from_past<-cbind(phylo$edge,node.age(phylo)$ages)
	nbtips<-Ntip(phylo)
	nbtips_tot<-nbtips_tot+nbtips
	ages<-rbind(from_past[,2:3],c(nbtips+1,0))
	ages<-ages[order(ages[,1]),]
	age<-max(ages[,2])
	#Here tj is the age of the ancestral split or the start times; sj1 is the stop time or the interval of 
	#daughter j1; and sj2 is the stop time or interval of daughter j2 until the next split. 
	#root == age.
	for (j in 1:(nbtips-1)){
		node<-(nbtips+j)
		edges<-phylo$edge[phylo$edge[,1]==node,]
		tj<-age-ages[edges[1,1],2]
		sj1<-age-ages[edges[1,2],2]
		sj2<-age-ages[edges[2,2],2]
		node.anc<-phy$edge[which(phy$edge[,2]==node),1]
		if (length(node.anc)==0) {
			node.anc<-node
		#assume params of root = params of root parent
		}
		else {
			ancestral.params<-rbind(ancestral.params, data.frame(node=node, 
																 turnover=getGeometricMeanParam(start.time=tj, stop.time=tj+phy$edge.length[j], param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$turnover, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k, param.splits=turnover.splits), 
																 eps=getGeometricMeanParam(start.time=tj, stop.time=tj+phy$edge.length[j], param.anc=ancestral.params[which(ancestral.params$node==node.anc),]$eps, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k, param.splits=eps.splits)))
		}
		indLikelihood<-c(indLikelihood,SetBirth(stop.time=tj, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)*Psi(s=sj1,t=tj, f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)*Psi(s=sj2,t=tj,f=f, turnover.param.anc=ancestral.params[which(ancestral.params$node==node),]$turnover, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=ancestral.params[which(ancestral.params$node==node),]$eps, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits))
	}
	#Can be removed -- we do not care about stem age
	#indLikelihood<-c(indLikelihood,Psi(s=age,t=tot_time,f=f,turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k))
	#Eq. 1 from Morlon et al 2011:
	data_lik<-prod(indLikelihood)*f^nbtips_tot
	Phi<-Phi(stop.time=tot_time,f=f, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	final_lik<-data_lik/(1-Phi)
	
	return(log(final_lik))
}

#sample likelihood
library(picante)
library(ape)
#phy<-rcoal(10)
phy<-read.tree("~/Dropbox/CollabBeaulieu/GeneralDiversification/test.tre")
#phy<-drop.tip(phy, paste("t", c(1:7), sep=""))
#phy<-drop.tip(phy, paste("t", c(1:5), sep=""))

#Rprof()
us.DD.LH<-(getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=0.5, turnover.sigma.time=0, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=1, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=0, eps.sigma.time=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, split.times=branching.times(phy), k=1000))
us.bd.LH<-(getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=0.5, turnover.sigma.time=0, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=0, eps.sigma.time=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, split.times=branching.times(phy), k=1000))
#Rprof(NULL)
#summaryRprof()

print(paste("us.bd.LH", us.bd.LH))
print(paste("us.DD.LH", us.DD.LH))
print(paste("us.DD.LH - us.bd.LH", us.DD.LH - us.bd.LH))
