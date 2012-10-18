##General diversification simulation##

library(ape)
source("SetParameter.R")

BirthTree<-function(phy, interval.length) {
	#do things: add a taxon, then run GrowTree
	lucky.taxon<-floor(runif(1,1,1+Ntip(phy)))
	phy<-bind.tree(x=phy, y=compute.brlen(stree(2,tip.label=c(lucky.taxon, Ntip(phy)+1)),method=0), where=lucky.taxon)
	phy<-GrowTree(phy, interval.length)
	return(phy)
}

DeathTree<-function(phy, interval.length) {
	#kill things: kill a taxon, then run GrowTree
	##drop.tip does not like pruning 1 from a 2 tip tree. Throws an error. This is my fix, not sure it is legal. 
	if(Ntip(phy)==2){
		return(phy=NULL)
	}
	else{
		unlucky.taxon<-floor(runif(1,1,1+Ntip(phy)))
		phy<-drop.tip(phy, unlucky.taxon)
		phy<-GrowTree(phy, interval.length)
	}
	return(phy)
}

GrowTree<-function(phy, interval.length) {
	#do things: lengthen terminal branches
	phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]+interval.length
	return(phy)
}

GetSim<-function(max.time=1, max.ntax=Inf, max.wall.time=Inf, check.file=NULL, start.file=NULL, prob.interval=0.0001, return.all.extinct=TRUE, verbose=TRUE, check.interval=1800, turnover.param.anc=1, turnover.sigma.indep=0, turnover.weight.anc=.53, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, eps.param.anc=.86, eps.sigma.indep=0, eps.weight.anc=1, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, k=Inf, warning.diversity=100) {
	#first initialize
	phy<-compute.brlen(stree(2, tip.label=c(2,1)),method=0)
	depth.time<-max.time
	start.time<-Sys.time()
	last.save.time<-Sys.time()
	split.times<-branching.times(phy)
	birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=split.times, eps.splits=split.times)
	death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=split.times, eps.splits=split.times)
	approx.expected.diversity=2*exp((birth-death)*max.time) #note that this is rough, as it does not take into account changing rates nor ascertainment bias
	if(approx.expected.diversity>warning.diversity) {
		warning(paste("Expected diversity is very roughly",approx.expected.diversity))
	}
	if(verbose) {
		print(paste("Expected diversity is very roughly",approx.expected.diversity))
	}
	interval.length<-qexp(prob.interval, rate=Ntip(phy)*(birth+death)) #how small should our interval be so that the chance of something happening in it is small
	print(paste("approx.expected.diversity*(birth+death)",approx.expected.diversity*(birth+death),approx.expected.diversity,(birth+death)))
	if(!is.null(start.file)) {
		load(start.file)
	}
	while(depth.time>0 & Ntip(phy)<=max.ntax & (Sys.time()-start.time)<max.wall.time) {
		depth.time<-depth.time-interval.length
		split.times<-branching.times(phy)
		birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=split.times, eps.splits=split.times)
		death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=split.times, eps.splits=split.times)
		if(rpois(1,Ntip(phy)*birth*interval.length)>0) {
			phy<-BirthTree(phy, interval.length)
			if(verbose) {
				print(c(Ntip(phy),depth.time))
			}
			interval.length<-qexp(prob.interval, rate=Ntip(phy)*(birth+death)) #how small should our interval be so that the chance of something happening in it is small
		}
		else if (rpois(1,Ntip(phy)*death*interval.length)>0) {
			phy<-DeathTree(phy, interval.length)
			#Moved this next bit up because then an error is thrown if tree returned from DeathTree is null:
			if (is.null(phy)) { 
				if(return.all.extinct) {
					return(NULL)
				}
				else {
					depth.time<-max.time
					phy<-compute.brlen(stree(2,tip.label=c(2,1)),method=0)
				}
			}
			interval.length<-qexp(prob.interval, rate=Ntip(phy)*(birth+death)) #how small should our interval be so that the chance of something happening in it is small
			if(verbose) {
				print(c(Ntip(phy),depth.time))
			}			
		}
		else {
			phy<-GrowTree(phy, interval.length)
		}
		if (Sys.time()-last.save.time>check.interval & !is.null(check.file)) {
			save(phy, depth.time, file=check.file)
			last.save.time<-Sys.time()
		}
	}
	if(verbose) {
		print(c(Ntip(phy),depth.time))
	}	
	phy<-GrowTree(phy, max.time-max(branching.times(phy)))
	phy$tip.label<-paste("t",c(1:Ntip(phy)),sep="")
	return(phy)
}

