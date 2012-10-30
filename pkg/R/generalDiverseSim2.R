##General diversification simulation##

library(ape)
source("SetParameter.R")

#our object while we are building up a tree is a data.frame
#first two columns are an edge matrix
#then edge length
#then boolean: T if  a tip (so it can keep growing
#then turnover in ancestor to that edge
#then turnover at end of that edge (just before the next split)
#then eps in ancestor to that edge
#then eps at end of that edge (just before the next split)
#as it grows, some taxa will go extinct. Those rows will be deleted
#at end, renumber so ape is happy, convert to phylo object
#need to pass turnover 

sim2phylo<-function(sim.object) {
	n.edge.all<-length(sim.object$to)
	n.tip<-sum(sim.object$tip)
	phy<-compute.brlen(stree(2, tip.label=c(2,1)),method=0)
	phy$edge<-matrix(NA,nrow=n.edge.all,ncol=2)
	phy$edge[,2]<-sequence(n.edge.all)
	equivalence<-data.frame(phylo=sequence(n.edge.all+1),sim.object=rep(NA,n.edge.all+1))
	equivalence[sequence(n.tip),2]<-sim.object$to[sim.object$tip]
	equivalence[c((n.tip+1):(n.edge.all+1)),2]<-sort(c(1,sim.object$to[!sim.object$tip]))
	print(equivalence)
	for(i in sequence(n.edge.all)) {
		relevant.row<-which(sim.object$to==equivalence[i,2])
		ancestor.node<-sim.object$from[relevant.row]
		phy$edge[i,1]<-equivalence[which(equivalence[,2]==ancestor.node), 1]
	}
	phy$edge.length<-sim.object$edge.length
	phy$Nnode<-1+n.edge.all-n.tip
	phy$tip.label=paste("t",sequence(n.tip),sep="")
	phy$order<-NULL
	phy<-reorder(phy)
	return(phy)
}

InitializeSimObject<-function(birth,death, turnover.anc, eps.anc) {
	return(data.frame(from=c(1,1),to=c(2,3), edge.length=rep(0,2), tip=rep(TRUE,2), turnover.anc=turnover.anc, turnover.present=rep(birth+death,2), eps.anc=eps.anc, eps.present=rep(death/birth,2)))
}

GrowSimObject<-function(sim.object, interval.length) {
	sim.object$edge.length[which(sim.object$tip)]<-sim.object$edge.length[which(sim.object$tip)] + interval.length
	return(sim.object)
}

UpdateTerminalParams<-function(sim.object) {
	#SetParameter for tips only
	Get eps.present, turnover.present, for all tips
}

BirthSimObject<-function(sim.object, interval.length, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, split.times, k, turnover.splits, turnover.sigma.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent,  eps.splits, eps.sigma.anc ) {
	lucky.tip<-(sim.object$to[which(sim.object$tip)])[floor(runif(1,1,1+sum(sim.object$tip)))]
	lucky.edge<-which(sim.object$to==lucky.tip)
	sim.object<-GrowSimObject(sim.object, interval.length)
	sim.object<-UpdateTerminalParams(sim.object)
	sim.object<-rbind(sim.object, sim.object[lucky.edge], sim.object[lucky.edge])
	sim.object$tip[lucky.edge]<-FALSE
	descendants<-c(dim(sim.object)[1]-1, dim(sim.object)[1])
	for (descendant.index in sequence(length(descendants))) {
		descendant<-descendants[descendant.index]
		sim.object$to[descendant]<-1+max(sim.object$to)
		sim.object$edge.length[descendant]<-0
		sim.object$turnover.anc[descendant] <- GetNewAncParam(stop.time, param.anc=sim.object$turnover.anc[lucky.edge], sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, param.splits=turnover.splits, param.sigma.anc=turnover.param.sigma.anc)
		sim.object$eps[descendant] <- GetNewAncParam(stop.time, param.anc=sim.object$eps.anc[lucky.edge], sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k, param.splits=eps.splits, param.sigma.anc=eps.param.sigma.anc)
	}
	
	
	sim.object<-UpdateTerminalParams(sim.object)
	return(sim.object)
}

DeathSimObject<-function(sim.object, interval.length) {
	unlucky.tip<-(sim.object$to[which(sim.object$tip)])[floor(runif(1,1,1+sum(sim.object$tip)))]
	unlucky.edge<-which(sim.object$to==unlucky.tip)
	sim.object<-GrowSimObject(sim.object, interval.length)
	merged.edge<-sim.object$from[which(sim.object$to==unlucky.tip)]
	surviving.descendant.of.merged.edge<-sim.object$to[which(sim.object$from==merged.edge)]
	surviving.descendant.of.merged.edge<-surviving.descendant.of.merged.edge[which(surviving.descendant.of.merged.edge!=unlucky.edge)]
	sim.object$from[surviving.descendant.of.merged.edge]<-sim.object$from[merged.edge]
	sim.object$edge.length[surviving.descendant.of.merged.edge]<-sim.object$edge.length[surviving.descendant.of.merged.edge]+sim.object$edge.length[merged.edge]
	sim.object<-sim.object[-unlucky.edge]
	sim.object<-UpdateTerminalParams(sim.object)
	return(sim.object)
}


BirthTree<-function(phy, interval.length) {
	#do things: add a taxon, then run GrowTree
	lucky.taxon<-floor(runif(1,1,1+Ntip(phy)))
	phy<-bind.tree(x=phy, y=compute.brlen(stree(2,tip.label=c(lucky.taxon, Ntip(phy)+1)),method=0), where=lucky.taxon)
	phy<-GrowTree(phy, interval.length)
	return(phy)
}

DeathTree<-function(phy, interval.length) {
	#kill things: kill a taxon, then run GrowTree
	##drop.tip does not like pruning 1 tip from a 2 tip tree. Throws an error. This is my fix for now. 
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

GetSim<-function(max.time=1, max.ntax=Inf, max.wall.time=Inf, check.file=NULL, start.file=NULL, prob.interval=0.001, return.all.extinct=TRUE, verbose=TRUE, check.interval=1800, turnover.param.indep=.7, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, eps.param.indep=.95, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, k=Inf, warning.diversity=100) {
	#first initialize
	#phy<-compute.brlen(stree(2, tip.label=c(2,1)),method=0)
	depth.time<-max.time
	start.time<-Sys.time()
	last.save.time<-Sys.time()
	
	#The split time = 0. Crashes.
	turnover.param.anc <- turnover.param.indep
	eps.param.anc <- eps.param.indep
	
	ancestral.params <- data.frame(node=Ntip(phylo)+1, turnover=turnover.param.anc, eps=eps.param.anc)
	turnover.splits <- exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep))
	eps.splits <- exp(rnorm(1, log(eps.param.indep), eps.sigma.indep))
	
	birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=c(0), k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=c(0), k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	approx.expected.diversity=2*exp((birth-death)*max.time) #note that this is rough, as it does not take into account changing rates nor ascertainment bias
	if(approx.expected.diversity>warning.diversity) {
		warning(paste("Expected diversity is very roughly",approx.expected.diversity))
	}
	if(verbose) {
		print(paste("Expected diversity is very roughly",approx.expected.diversity))
	}
	sim.object<-InitializeSimObject(birth,death, turnover.param,anc, eps.param.anc)
	phy<-sim2phylo(sim.object)
	split.times<-branching.times(phy)+depth.time

	interval.length<-qexp(prob.interval, rate=Ntip(phy)*(birth+death)) #how small should our interval be so that the chance of something happening in it is small
	print(paste("approx.expected.diversity*(birth+death)",approx.expected.diversity*(birth+death),approx.expected.diversity,(birth+death)))
	if(!is.null(start.file)) {
		load(start.file)
	}
	while(depth.time>0 & Ntip(phy)<=max.ntax & (Sys.time()-start.time)<max.wall.time) {
#		depth.time<-depth.time-interval.length #change to depth.time<-depth.time-rexp(1, 1/turnover rate (check, might be turnover rate))
		interval.length<-rexp(1, 1/(birth+death))
		depth.time<-depth.time-interval.length
		if (depth.time<0) {
			phy<-GrowTree(phy, max.time-max(branching.times(phy)))
			phy$tip.label<-paste("t",c(1:Ntip(phy)),sep="")
			return(phy)
		}
		split.times<-branching.times(phy)+depth.time
		position <- max(which(split.times>=depth.time),1)
		stop.time <- mean(split.times[position], depth.time)
				
		birth<-SetBirth(stop.time=max.time-depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=turnover.splits, eps.splits=eps.splits)
		death<-SetDeath(stop.time=max.time-depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=turnover.splits, eps.splits=eps.splits)

		if(runif(1,0,1)<(birth/(birth+death))) {
			turnover.splits <- c(rep(NA,length(turnover.splits)),exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep)))
			eps.splits <- c(rep(NA,length(eps.splits)),exp(rnorm(1, log(eps.param.indep), eps.sigma.indep)))

			#Change to BirthSimObject here:
			phy<-BirthTree(phy, interval.length)
			#
			#Then call remake phy:
			#
			if(verbose) {
				print(c(Ntip(phy),depth.time))
			}
		}
		else {
			#Change to DeathSimObject here:
			#Describe:
			turnover.splits<-turnover.splits[-1]
			eps.splits<-eps.splits[-1]
			phy<-DeathTree(phy, interval.length)
			#
			#Then call remake phy:
			#
			if (is.null(phy)) { 
				if(return.all.extinct) {
					return(NULL)
				}
				else {
					depth.time<-max.time
					#Just reinitialize here:
					phy<-compute.brlen(stree(2,tip.label=c(2,1)),method=0)
					#
				}
			}
			if(verbose) {
				print(c(Ntip(phy),depth.time))
			}			
		}
		if (Sys.time()-last.save.time>check.interval & !is.null(check.file)) {
			save(phy, depth.time, file=check.file)
			last.save.time<-Sys.time()
		}
	}
}

