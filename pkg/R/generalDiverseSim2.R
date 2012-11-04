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

GetSim<-function(max.time=1, max.ntax=Inf, max.wall.time=Inf, check.file=NULL, start.file=NULL, prob.interval=0.001, return.all.extinct=TRUE, verbose=TRUE, check.interval=1800, turnover.param.indep=.55, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=1, turnover.trend.scaling=0, turnover.trend.exponent=0, eps.param.indep=0.85, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, k=40, turnover.sigma.anc=0, eps.sigma.anc=0, warning.diversity=100) {

	depth.time<-max.time
	start.time<-Sys.time()
	last.save.time<-Sys.time()
	options(digits=10)
	
	turnover.param.anc <- turnover.param.indep
	eps.param.anc <- eps.param.indep
	
	turnover.splits <- exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep))
	eps.splits <- exp(rnorm(1, log(eps.param.indep), eps.sigma.indep))
	
	birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=c(0), k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=c(0), k, turnover.splits=turnover.splits, eps.splits=eps.splits)

	approx.expected.diversity=2*exp((birth-death)*max.time) #note that this is rough, as it does not take into account changing rates nor ascertainment bias
	
	if(approx.expected.diversity>warning.diversity) {
		warning(paste("Expected diversity is very roughly", approx.expected.diversity))
	}
	if(verbose) {
		print(paste("Expected diversity is very roughly", approx.expected.diversity))
	}
	
	sim.object<-InitializeSimObject(birth,death, turnover.param.anc, eps.param.anc)
	phy<-sim2phylo(sim.object)
	split.times<-branching.times(phy)+depth.time
	
	print(paste("approx.expected.diversity*(birth+death)",approx.expected.diversity*(birth+death),approx.expected.diversity,(birth+death)))
	
	if(!is.null(start.file)) {
		load(start.file)
	}
	
	while(depth.time>0 & Ntip(phy)<=max.ntax & (Sys.time()-start.time)<max.wall.time) {

		alive<-length(branching.times(phy))
		interval.length<-rexp(1, alive*(birth+death))
		depth.time<-depth.time-interval.length

		if (depth.time<0) {
			
			###########################REMAINING ISSUE###########################
			#If node 1 remains in sim.object after max.time, then the total height should equal max.time, right?
			phy <- sim2phylo(sim.object)
			sim.object <- GrowSimObject(sim.object, max.time-max(branching.times(phy)))
			#####################################################################
			
			phy <- sim2phylo(sim.object)

			###########################REMAINING ISSUE###########################
			phy <- reorder(phy,"pruningwise")
			phy <- reorder(phy,"cladewise")
			#####################################################################

			return(phy)
		}
		
		split.times<-sort(branching.times(phy)+depth.time, decreasing=TRUE)
		
		###########################REMAINING ISSUE###########################
		#What we had before -- cannot remember why we use NAs:
#		turnover.splits <- c(rep(NA,length(turnover.splits)),exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep)))
#		eps.splits <- c(rep(NA,length(eps.splits)),exp(rnorm(1, log(eps.param.indep), eps.sigma.indep)))
		#The following should be equivalent. The benefit of doing it this way, is that we do not have to worry about position:
###Should not be NA:
		turnover.splits <- rep(exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep)),length(split.times))
		eps.splits <- rep(exp(rnorm(1, log(eps.param.indep), turnover.sigma.indep)),length(split.times))
		#####################################################################

		birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=turnover.splits, eps.splits=eps.splits)
		death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, split.times=split.times, k, turnover.splits=turnover.splits, eps.splits=eps.splits)
		
		birth.proportion<-0
		if((birth+death)>0) {
			birth.proportion<-birth/(birth+death)
		}
		if(runif(1,0,1)<birth.proportion){
			sim.object <- BirthSimObject(sim.object=sim.object, interval.length=interval.length, stop.time=depth.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
			phy <- sim2phylo(sim.object)
			if(verbose) {
				print(c(Ntip(phy),depth.time))
			}
		}
		else {
			sim.object <- DeathSimObject(sim.object=sim.object, interval.length=interval.length, stop.time=depth.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)

			###########################REMAINING ISSUE###########################
			#Necessary? Probably not, because again, the vector is being remade each time
			#but is that right?
			turnover.splits<-turnover.splits[-1]
			eps.splits<-eps.splits[-1]
			#####################################################################

			if(!is.null(sim.object)){
				phy <- sim2phylo(sim.object) 
			}
			if (is.null(sim.object)) { 
				if(return.all.extinct) {
					return(NULL)
				}
				else {
					depth.time<-max.time
					sim.object<-InitializeSimObject(birth, death, turnover.param.indep, eps.param.indep)
					phy <- sim2phylo(sim.object)
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

sim2phylo<-function(sim.object) {
	n.edge.all<-length(sim.object$to)
	n.tip<-sum(sim.object$tip)
	phy<-compute.brlen(stree(2, tip.label=c(2,1)),method=0)
	phy$edge<-matrix(NA,nrow=n.edge.all,ncol=2)
	equivalence<-data.frame(phylo=sequence(n.edge.all+1),sim.object=rep(NA,n.edge.all+1))
	equivalence[sequence(n.tip),2]<-sim.object$to[sim.object$tip]
	#sim.object$from[which(!sim.object$from%in%sim.object$to)][1] is to designate the root in case "1" is murdered:
	print(sim.object$from[which(!sim.object$from%in%sim.object$to)])
	print(sim.object)
	equivalence[c((n.tip+1):(n.edge.all+1)),2]<-sort(c(sim.object$from[which(!sim.object$from%in%sim.object$to)][1],sim.object$to[!sim.object$tip]))
	for(i in sequence(n.edge.all)) {
		#Defines the row in the equivalence table that is equal to the row in the sim.object:
		phy$edge[i,2]<-equivalence[which(sim.object$to[i]==equivalence[,2]),1]
		phy$edge[i,1]<-equivalence[which(sim.object$from[i]==equivalence[,2]),1]
	}
	phy$edge.length<-sim.object$edge.length
	phy$Nnode<-1+n.edge.all-n.tip
	phy$tip.label=paste("t",sequence(n.tip),sep="")
	phy$order<-NULL
	phy<-reorder(phy)
	return(phy)
}

InitializeSimObject<-function(birth, death, turnover.anc, eps.anc) {
	return(data.frame(from=c(1,1),to=c(2,3), edge.length=rep(0,2), tip=rep(TRUE,2), turnover.anc=turnover.anc, turnover.present=rep(birth+death,2), eps.anc=eps.anc, eps.present=rep(death/birth,2)))
}

GrowSimObject<-function(sim.object, interval.length) {
	sim.object$edge.length[which(sim.object$tip)]<-sim.object$edge.length[which(sim.object$tip)] + interval.length
	return(sim.object)
}

UpdateTerminalParams<-function(sim.object, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, split.times, k, turnover.splits, turnover.sigma.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent,  eps.splits, eps.sigma.anc) {
	descendants<-c(dim(sim.object)[1]-1, dim(sim.object)[1])
	for (descendant.index in sequence(length(descendants))) {
		descendant<-descendants[descendant.index]
		if(sim.object$tip[descendant]==TRUE){
			sim.object$turnover.present[descendant] <- SetParameter(stop.time, param.anc=sim.object$turnover.anc[descendant], sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, param.splits=turnover.splits)
			sim.object$eps.present[descendant] <- SetParameter(stop.time, param.anc=sim.object$eps.anc[descendant], sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k, param.splits=eps.splits)
		}
	}
	return(sim.object)
}

BirthSimObject<-function(sim.object, interval.length, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, split.times, k, turnover.splits, turnover.sigma.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, eps.splits, eps.sigma.anc){
	lucky.tip<-(sim.object$to[which(sim.object$tip)])[floor(runif(1,1,1+sum(sim.object$tip)))]
	lucky.edge<-which(sim.object$to==lucky.tip)
	sim.object<-GrowSimObject(sim.object, interval.length)
	sim.object<-UpdateTerminalParams(sim.object=sim.object, stop.time=stop.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
	sim.object<-rbind(sim.object, sim.object[lucky.edge,], sim.object[lucky.edge,])
	sim.object$tip[lucky.edge]<-FALSE
	descendants<-c(dim(sim.object)[1]-1, dim(sim.object)[1])
	for (descendant.index in sequence(length(descendants))) {
		descendant<-descendants[descendant.index]
		sim.object$to[descendant]<-1+max(sim.object$to)
		sim.object$from[descendant]<-lucky.tip
		sim.object$edge.length[descendant]<-0
		sim.object$turnover.anc[descendant] <- GetNewAncParam(stop.time, param.anc=sim.object$turnover.anc[lucky.edge], sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.scaling=turnover.trend.scaling, trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, param.splits=turnover.splits, param.sigma.anc=turnover.sigma.anc)
		sim.object$eps.anc[descendant] <- GetNewAncParam(stop.time, param.anc=sim.object$eps.anc[lucky.edge], sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.scaling=eps.trend.scaling, trend.exponent=eps.trend.exponent, split.times=split.times, k=k, param.splits=eps.splits, param.sigma.anc=eps.sigma.anc)
	}
	sim.object<-UpdateTerminalParams(sim.object=sim.object, stop.time=stop.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
	return(sim.object)
}

DeathSimObject<-function(sim.object, interval.length, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.scaling, turnover.trend.exponent, split.times, k, turnover.splits, turnover.sigma.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.scaling, eps.trend.exponent, eps.splits, eps.sigma.anc) {
	
	###########################REMAINING ISSUE###########################
	#Error is thrown when N<=2 during a death event -- not ideal. Although, 
	#we would have to live with this if we were still relying on ape...
	if(length(sim.object[,1])<=2){
		return(NULL)
	}
	#####################################################################
	else{
		unlucky.tip<-(sim.object$to[which(sim.object$tip)])[floor(runif(1,1,1+sum(sim.object$tip)))]
		sim.object<-GrowSimObject(sim.object, interval.length)
		#Defines the new tip:
		merged.edge<-sim.object$from[which(sim.object$to==unlucky.tip)]
		#Defines the unlucky tip and its sister:
		surviving.descendant.of.merged.edge<-sim.object$to[which(sim.object$from==merged.edge)]
		#Defines the sister of the tip:
		surviving.descendant.of.merged.edge<-surviving.descendant.of.merged.edge[which(surviving.descendant.of.merged.edge!=unlucky.tip)]
		#Combines the pieces of the two merged edges together if the sister is also a tip:
		if(sim.object$tip[which(sim.object$to==surviving.descendant.of.merged.edge)]==TRUE){			
			sim.object$edge.length[which(sim.object$to==merged.edge)]<-sim.object$edge.length[which(sim.object$to==surviving.descendant.of.merged.edge)]+sim.object$edge.length[which(sim.object$to==merged.edge)]
			#Defines the ancestor of the unlucky tip
			unlucky.edge<-sim.object$to[which(sim.object$from==merged.edge)]
			#Drops the necessary rows:
			sim.object<-sim.object[-which(sim.object$to==unlucky.edge),]
			#Changes the new tip to a tip:
			sim.object$tip[which(sim.object$to==merged.edge)]<-TRUE
		}
		#If the sister is not a tip then the to must be changed for the ancestor of the non-tip sister
		else{
			unlucky.edge<-sim.object$to[which(sim.object$from==merged.edge)]
			to.be.deleted<-which(sim.object$to==unlucky.edge)
			sim.object$to[which(sim.object$to==merged.edge)]<-surviving.descendant.of.merged.edge
			new.ancestor<-sim.object$from[which(sim.object$to==surviving.descendant.of.merged.edge)]
			sim.object$from[which(sim.object$to==surviving.descendant.of.merged.edge)]<-min(new.ancestor)
			sim.object$edge.length[which(sim.object$to==surviving.descendant.of.merged.edge)]<-sum(sim.object$edge.length[which(sim.object$to==surviving.descendant.of.merged.edge)])
			sim.object<-sim.object[-to.be.deleted,]
		}
		#Calculate the current rates at all the tips:
		sim.object<-UpdateTerminalParams(sim.object=sim.object, stop.time=stop.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.scaling=turnover.trend.scaling, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, k=k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.scaling=eps.trend.scaling, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
		return(sim.object)
	}
}

#sim.object2<-BirthSimObject(sim.object, interval.length=.5, stop.time=8, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, split.times=branching.times(phy), k=1, turnover.splits=c(rep(.7,length(branching.times(phy)))), turnover.sigma.anc=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, eps.splits=c(rep(.95,length(branching.times(phy)))), eps.sigma.anc=0)

#sim.object2<-BirthSimObject(sim.object2, interval.length=.5, stop.time=7, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, split.times=branching.times(phy), k=1, turnover.splits=c(rep(.7,length(branching.times(phy)))), turnover.sigma.anc=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, eps.splits=c(rep(.95,length(branching.times(phy)))), eps.sigma.anc=0)

#sim.object2<-BirthSimObject(sim.object2, interval.length=.5, stop.time=6, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, split.times=branching.times(phy), k=1, turnover.splits=c(rep(.7,length(branching.times(phy)))), turnover.sigma.anc=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, eps.splits=c(rep(.95,length(branching.times(phy)))), eps.sigma.anc=0)

#sim.object2<-BirthSimObject(sim.object2, interval.length=.5, stop.time=5.5, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, split.times=branching.times(phy), k=1, turnover.splits=c(rep(.7,length(branching.times(phy)))), turnover.sigma.anc=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, eps.splits=c(rep(.95,length(branching.times(phy)))), eps.sigma.anc=0)

#sim.object3<-DeathSimObject(sim.object2, interval.length=.5, stop.time=5, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, split.times=branching.times(phy), k=1, turnover.splits=c(rep(.7,length(branching.times(phy)))), turnover.sigma.anc=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, eps.splits=c(rep(.95,length(branching.times(phy)))), eps.sigma.anc=0)

#sim.object4<-DeathSimObject(sim.object3, interval.length=.5, stop.time=5, turnover.sigma.indep=0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.scaling=0, turnover.trend.exponent=0, split.times=branching.times(phy), k=1, turnover.splits=c(rep(.7,length(branching.times(phy)))), turnover.sigma.anc=0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.scaling=0, eps.trend.exponent=0, eps.splits=c(rep(.95,length(branching.times(phy)))), eps.sigma.anc=0)

