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
#	split.times = the vector of split times
#	k = carrying capacity

#general.diversity<-function(, params.turnover.constant, params.eps.constant){

#Inputs to getLikelihood...:(phylo,tot_time,f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k)

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

SetLogParameter <- function(current.time, param.indep, param.anc, sigma.time, sigma.indep, weight.anc, weight.logistic, trend.scaling=0, trend.exponent=1, split.times, k) {
	
	if (trend.exponent<0) {
		warning("It is possible to have a trend.exponent less than zero, but this makes little sense")
	}
	n.taxa<-1+max(which(split.times>=current.time),1)
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc)
	param.mean <- (param.starting + trend.scaling * (current.time ^ trend.exponent)) * (1 - (weight.logistic * n.taxa / k)) 
	return(list(expected = param.mean, random = rnorm(1, param.mean, sigma.time * current.time + sigma.indep)))
} 

SetBirth <- function(current.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	turnover<-lapply(SetLogParameter(current.time=current.time, param.indep=turnover.param.indep, param.anc=turnover.param.anc, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k), exp)
	extinction.fraction<-lapply(SetLogParameter(current.time=current.time, param.indep=eps.param.indep, param.anc=eps.param.anc, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k), exp)
	birth.expected<-turnover$expected / (1 + extinction.fraction$expected)
	birth.random<-turnover$random / (1 + extinction.fraction$random)
	return(list(expected=birth.expected, random=birth.random))
}

SetDeath <- function(current.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	turnover<-lapply(SetLogParameter(current.time=current.time, param.indep=turnover.param.indep, param.anc=turnover.param.anc, sigma.time=turnover.sigma.time, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k), exp)
	extinction.fraction<-lapply(SetLogParameter(current.time=current.time, param.indep=eps.param.indep, param.anc=eps.param.anc, sigma.time=eps.sigma.time, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k), exp)
	death.expected<-(extinction.fraction$expected * turnover$expected) / (1 + extinction.fraction$expected)
	death.random<-(extinction.fraction$random * turnover$random) / (1 + extinction.fraction$random)
	return(list(expected=death.expected, random=death.random))
}

SetDiversification <- function(current.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	return(list(expected=(SetBirth(current.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$expected - SetDeath(current.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$expected), 
				random=(SetBirth(current.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$random - SetDeath(current.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$random)))
}

SetDiversificationExpected <- function(current.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	return(SetDiversification(current.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$expected)
}

IntegrateDiversificationOverTime <- function(stop.time, start.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	return(integrate(Vectorize(SetDiversificationExpected, vectorize.args="current.time"), lower=stop.time, upper=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, stop.on.error=FALSE)$value)
}

IntegrateDiversificationOverTime.int.0 <- function(start.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	result<-IntegrateDiversificationOverTime(stop.time=0, start.time=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)
	lambda<-SetBirth(current.time=start.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$expected
	result<-exp(result)*lambda
	return(result)
}

IntegrateDiversificationOverTime.int.int <- function(stop.time, start.time, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	result<-(integrate(Vectorize(IntegrateDiversificationOverTime.int.0, vectorize.args="start.time"), lower=stop.time, upper=start.time,  turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k, stop.on.error=FALSE)$value)
	return(result)
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

Phi<-function(current.time, f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	res<-(1-exp(IntegrateDiversificationOverTime(stop.time=0,start.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)) / (1/f+IntegrateDiversificationOverTime.int.int(stop.time=0,start.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)))
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

Psi<-function(s, current.time, f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k) {
	#Equation broken up for ease of debugging. The extreme number of input makes it difficult to see what is going on:
	#res = first part of Morlon et al 2011 eq. 3
	res<-exp(IntegrateDiversificationOverTime(stop.time=s, start.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k))
	#res2 = numerator of Morlon et al 2011 eq. 3
	res2<-IntegrateDiversificationOverTime.int.int(stop.time=s, start.time=current.time, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)
	#res3 = inside brackets of Morlon et al 2011 eq. 3
	res3<-res2/(1/f+IntegrateDiversificationOverTime.int.int(stop.time=0, start.time=s, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k))
	#res = Morlon et al 2011 eq. 3
	res<-res*(abs(1+res3)^(-2))
	return(res)
}

######################################################################################################################################
######################################################################################################################################
### Obtains the likelihood:
######################################################################################################################################
######################################################################################################################################

getLikelihood.gen.bd<-function(phylo,tot_time,f, turnover.param.indep, turnover.param.anc, turnover.sigma.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.indep, eps.param.anc, eps.sigma.time, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times, k){
	
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
		indLikelihood<-c(indLikelihood,SetBirth(current.time=tj, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)$expected*Psi(s=sj1,current.time=tj,f=f, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)*Psi(s=sj2,current.time=tj,f=f, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k))
		#Takes about 60 seconds per node#
		#print("here")
	}
	indLikelihood<-c(indLikelihood,Psi(s=age,current.time=tot_time,f=f,turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k))
	#Eq. 1 from Morlon et al 2011:
	data_lik<-prod(indLikelihood)*f^nbtips_tot
	Phi<-Phi(current.time=tot_time,f=f, turnover.param.indep=turnover.param.indep, turnover.param.anc=turnover.param.anc, turnover.sigma.time=turnover.sigma.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, eps.param.indep=eps.param.indep, eps.param.anc=eps.param.anc, eps.sigma.time=eps.sigma.time, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, split.times=split.times, k=k)
	final_lik<-data_lik/(1-Phi)
	
	return(log(final_lik))
}


#sample likelihood
library(picante)
library(ape)
#phy<-rcoal(10)
phy<-read.tree("~/Dropbox/CollabBeaulieu/GeneralDiversification/test.tre")

print(getLikelihood.gen.bd(phylo=phy,tot_time=max(branching.times(phy)),f=1, turnover.param.indep=log(.75), turnover.param.anc=log(.75), turnover.sigma.time=.1, turnover.sigma.indep=.1, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=1, eps.param.indep=log(0.5), eps.param.anc=log(0.5), eps.sigma.time=.1, eps.sigma.indep=.1, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=1, split.times=branching.times(phy), k=100))
