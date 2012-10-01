#Overall goal is to have turnover rate or extinction fraction vary as a result of various
#   parameters.
#Note that with logistic growth, it makes the most sense to pass in n.taxa as actual
#   n.taxa (including those with no descendants, but in practice, n.taxa reconstructed
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
#	n.taxa = the vector of split times
#	k = carrying capacity

#general.diversity<-function(, params.turnover.constant, params.eps.constant){

#Inputs to getLikelihood...:(phylo,tot_time,f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k)

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

SetLogParameter <- function(stop.time, param.indep, param.anc, sigma.time, sigma.indep, weight.anc, weight.logistic, trend.scaling=0, trend.exponent=1, n.taxa, k) {
	
	if (trend.exponent<0) {
		warning("It is possible to have a trend.exponent less than zero, but this makes little sense")
	}
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc)
	param.mean <- (param.starting + trend.scaling * (stop.time ^ trend.exponent)) * (1 - (weight.logistic * n.taxa / k)) 
	return(list(expected = param.mean, random = exp(rnorm(1, log(param.mean), sigma.time * stop.time + sigma.indep))))
} 

SetBirth <- function(stop.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	turnover<-SetLogParameter(stop.time=stop.time, param.indep=turnover.param.indep, param.anc=turnover.param.anc, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, n.taxa=n.taxa, k=k)
	extinction.fraction<-SetLogParameter(stop.time=stop.time, param.indep=eps.param.indep, param.anc=eps.param.anc, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)
	birth.expected<-turnover$expected / (1 + extinction.fraction$expected)
	birth.random<-turnover$random / (1 + extinction.fraction$random)
	return(list(expected=birth.expected, random=birth.random))
}

SetDeath <- function(stop.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	turnover<-SetLogParameter(stop.time=stop.time, param.indep=turnover.param.indep, param.anc=turnover.param.anc, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, n.taxa=n.taxa, k=k)
	extinction.fraction<-SetLogParameter(stop.time=stop.time, param.indep=eps.param.indep, param.anc=eps.param.anc, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)
	death.expected<-(extinction.fraction$expected * turnover$expected) / (1 + extinction.fraction$expected)
	death.random<-(extinction.fraction$random * turnover$random) / (1 + extinction.fraction$random)
	return(list(expected=death.expected, random=death.random))
}

SetDiversification <- function(stop.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	return(list(expected=(SetBirth(stop.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$expected - SetDeath(stop.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$expected), 
				random=(SetBirth(stop.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$random - SetDeath(stop.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$random)))
}

SetDiversificationExpected <- function(stop.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	diversification<-SetDiversification(stop.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$expected
	return(diversification)
}

#Not needed right now, but going to keep in case the BM step requires this help:
#IntegrateSetDiversificationExpected <- function(lower, upper, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
#	if (turnover.weight.logistic>0 | eps.weight.logistic>0) {
#		boundaries<-n.taxa[which(n.taxa>lower)]
#		boundaries<-boundaries[which(boundaries<upper)]
#		lower.bounds<-c(lower, boundaries)
#		upper.bounds<-c(boundaries, upper)
#		integration.result<-0
#		for (i in sequence(length(lower.bounds))) {
#			integration.result<-integration.result + integrate(Vectorize(SetDiversificationExpected), lower=lower.bounds[i], upper=upper.bounds[i], turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k, stop.on.error=FALSE)$value
#		}
#		return(integration.result)
#	}
#	if (brownian motion) {
	
#	}
#	else {
#	  return(integrate(Vectorize(SetDiversificationExpected), lower=lower, upper=upper, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k, stop.on.error=FALSE)$value)
#	}
#}

IntegrateDiversificationOverTime <- function(stop.time, start.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	return(integrate(Vectorize(SetDiversificationExpected), lower=stop.time, upper=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k, stop.on.error=FALSE)$value)
#	return(IntegrateSetDiversificationExpected(lower=stop.time, upper=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k))
}

##The arguments here do not match below, obviously. But in Morlon, not a problem. Here a problem. And setting the vectorize.arg SLLLOOOWWWS things down.
IntegrateDiversificationOverTime.int.0 <- function(start.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	return(exp(IntegrateDiversificationOverTime(stop.time=0, start.time=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k))*SetBirth(stop.time=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$expected)
}

IntegrateDiversificationOverTime.int.int <- function(stop.time, start.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	return(integrate(Vectorize(IntegrateDiversificationOverTime.int.0), lower=stop.time, upper=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k, stop.on.error=FALSE)$value)
}

######################################################################################################################################
######################################################################################################################################
### This function computes Phi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

#Phi<-function(t,f.lamb,f.mu,f){
#	r<-function(x){f.lamb(x)-f.mu(x)}	
#	r.int<-function(x,y){integrate(Vectorize(r),x,y,stop.on.error=FALSE)$value}
#	r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
#	r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
#	res<-1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
#	res<-exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)

#	return(res)
#}

Phi<-function(stop.time, f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	res<-(1-exp(IntegrateDiversificationOverTime(stop.time=0,start.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)) / (1/f+IntegrateDiversificationOverTime.int.int(stop.time=0,start.time=stop.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)))
	return(res)	
}

######################################################################################################################################
######################################################################################################################################
### This function computes Psi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

#Psi<-function(s,t,f.lamb,f.mu,f, params.turnover.constant, params.eps.constant, ancestral.t){
#	r<-function(x){f.lamb(x)-f.mu(x)}	
#	r.int<-function(x,y){integrate(Vectorize(r),x,y,stop.on.error=FALSE)$value}
#	r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
#	r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}

#	res<-exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
#	return(res)
#}

Psi<-function(s, t, f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k) {
	#Equation broken up for ease of debugging. The extreme number of input makes it difficult to see what is going on:
	#res = first part of Morlon et al 2011 eq. 3
	res<-exp(IntegrateDiversificationOverTime(stop.time=s, start.time=t, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k))
	#res2 = numerator of Morlon et al 2011 eq. 3
	res2<-IntegrateDiversificationOverTime.int.int(stop.time=s, start.time=t, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)
	#res3 = inside brackets of Morlon et al 2011 eq. 3
	res3<-res2/(1/f+IntegrateDiversificationOverTime.int.int(stop.time=0, start.time=s, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k))
	#res = Morlon et al 2011 eq. 3
	res<-res*(abs(1+res3)^(-2))
	return(res)
}

######################################################################################################################################
######################################################################################################################################
### Obtains the likelihood:
######################################################################################################################################
######################################################################################################################################

getLikelihood.gen.bd<-function(phylo,tot_time,f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, n.taxa, k){
	
	indLikelihood<-c()
	nbfinal<-length(phylo)
	nbtips_tot<-0
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
		#Added n.taxa here, because when it is set above, time is no longer the same, and therefore the diversity is always equal to 2. I think this is what
		#was messing up the integration and causing it slow down.
		n.taxa<-1+max(which(branching.times(phylo)>=tj),1)
		indLikelihood<-c(indLikelihood,SetBirth(stop.time=tj, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)$expected*Psi(s=sj1,t=tj, f=f, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)*Psi(s=sj2,t=tj,f=f, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k))
	}
	indLikelihood<-c(indLikelihood,Psi(s=age,t=tot_time,f=f,turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k))
	#Eq. 1 from Morlon et al 2011:
	data_lik<-prod(indLikelihood)*f^nbtips_tot
	Phi<-Phi(stop.time=tot_time,f=f, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, n.taxa=n.taxa, k=k)
	final_lik<-data_lik/(1-Phi)
	
	return(log(final_lik))
}


#sample likelihood
library(picante)
library(ape)
#phy<-rcoal(10)
phy<-read.tree("~/Dropbox/CollabBeaulieu/GeneralDiversification/test.tre")
phy<-drop.tip(phy, paste("t", c(1:7), sep=""))

us.DD.LH<-(getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=.5, turnover.param.anc=.5, turnover.sigma.time=.1, turnover.sigma.indep=.1, turnover.weight.anc=0, turnover.weight.logistic=1, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=0, eps.param.anc=0, eps.sigma.time=.1, eps.sigma.indep=.1, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, n.taxa=branching.times(phy), k=100))

us.bd.LH<-(getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=.5, turnover.param.anc=.5, turnover.sigma.time=.1, turnover.sigma.indep=.1, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=0, eps.param.anc=0, eps.sigma.time=.1, eps.sigma.indep=.1, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, n.taxa=branching.times(phy), k=100))

print(paste("us.bd.LH", us.bd.LH))
print(paste("us.DD.LH", us.DD.LH))
print(paste("us.DD.LH - us.bd.LH", us.DD.LH - us.bd.LH))
