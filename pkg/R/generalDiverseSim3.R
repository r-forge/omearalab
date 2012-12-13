##General diversification simulation##

library(ape)
#source("SetParameter.R")

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

GetSim<-function(max.time=1, max.ntax=Inf, max.wall.time=Inf, check.file=NULL, start.file=NULL, return.all.extinct=TRUE, verbose=FALSE, check.interval=1800, turnover.param.indep=.10, turnover.sigma.indep=0.0, turnover.weight.anc=0, turnover.weight.logistic=0, turnover.trend.exponent=0, turn.k=Inf, turnover.sigma.anc=0, eps.param.indep=0.0, eps.sigma.indep=0, eps.weight.anc=0, eps.weight.logistic=0, eps.trend.exponent=0, eps.k=Inf, eps.sigma.anc=0, warning.diversity=Inf) {
	
	depth.time<-max.time
	start.time<-Sys.time()
	last.save.time<-Sys.time()
		options(digits=10)
	root.tracker=matrix(c(1,depth.time),nrow=1,ncol=2)
	turnover.param.anc <- turnover.param.indep
	eps.param.anc <- eps.param.indep
	
	turnover.splits <- exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep))
	eps.splits <- exp(rnorm(1, log(eps.param.indep), eps.sigma.indep))
	
	rate.track<-turnover.splits[1]
	time.track<-depth.time
	
	birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times=c(0), turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times=c(0), turn.k, eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
	
	#note that this is rough, as it does not take into account changing rates nor ascertainment bias
	approx.expected.diversity=2*exp((birth-death)*max.time) 	
	
	if(is.finite(max.time)){
		if(approx.expected.diversity>warning.diversity) {
			warning(paste("Expected diversity is very roughly", approx.expected.diversity))
		}
		if(verbose) {
			print(paste("Expected diversity is very roughly", approx.expected.diversity))
		}
	}
	
	sim.object<-InitializeSimObject(birth,death, turnover.param.anc, eps.param.anc)
	max.tip.count<-max(sim.object$to)
	phy<-sim2phylo(sim.object)
	split.times<-branching.times(phy)+depth.time
	
	if(is.finite(max.time)){
		print(paste("approx.expected.diversity*(birth+death)",approx.expected.diversity*(birth+death),approx.expected.diversity,(birth+death)))
	}
	
	if(!is.null(start.file)) {
		load(start.file)
	}
	
	while(depth.time>0 & Ntip(phy)<=max.ntax & (Sys.time()-start.time)<max.wall.time) {
		
		alive<-Ntip(phy)
		
		##Should be the sum of all rates (b and d) -- the ntax is now unnecessary:
#		interval.length<-rexp(1, alive*(birth+death))
		interval.length<-rexp(1, sum(sim.object$turnover.anc[which(sim.object$tip==TRUE)],sim.object$eps.anc[which(sim.object$tip==TRUE)]))
		
		depth.time<-depth.time-interval.length
		
		if(is.finite(max.time)){
			if (depth.time<0) {
				root.node<-sim.object$from[which(!sim.object$from%in%sim.object$to)][1]
				root.depth<-root.tracker[which(root.tracker[,1]==root.node),2]
				phy <- sim2phylo(sim.object)
				sim.object <- GrowSimObject(sim.object, root.depth-max(branching.times(phy)))
				print(sim.object)
				phy <- sim2phylo(sim.object)
				###########################REMAINING ISSUE###########################
				phy <- reorder(phy,"pruningwise")
				phy <- reorder(phy,"cladewise")
				#####################################################################
				return(phy)
			}
		}
		else{
			phy <- sim2phylo(sim.object)
			if(Ntip(phy)==max.ntax){
				obj<-NULL
				###########################REMAINING ISSUE###########################
				phy <- reorder(phy,"pruningwise")
				phy <- reorder(phy,"cladewise")
				#####################################################################
#				obj$phy<-phy
#				obj$rates<-rate.track
				obj$sim.object<-sim.object
				return(obj)
#				return(phy)
			}
		}
		split.times<-sort(branching.times(phy)+depth.time, decreasing=TRUE)
		turnover.splits <- rep(exp(rnorm(1, log(turnover.param.indep), turnover.sigma.indep)), length(split.times)) #gets refilled each interval, even though the last one is the only one used
		eps.splits <- rep(exp(rnorm(1, log(eps.param.indep), eps.sigma.indep)), length(split.times))
		
		birth<-SetBirth(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
		death<-SetDeath(stop.time=depth.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times=split.times, turn.k, eps.k, turnover.splits=turnover.splits, eps.splits=eps.splits)
				
		the.chosen.one<-which(rmultinom(1,1,prob=c(sim.object$turnover.anc[which(sim.object$tip==TRUE)],sim.object$eps.anc[which(sim.object$tip==TRUE)]))==1)
		#Choose a tip using a rmultinom with the first half of the prob vector corresponding to birth, and the second half corresponding to death.
		#If the lucky.tip is less than or equal to the number of tips, then it is a birth event:
		if(the.chosen.one <= alive){
#		if(runif(1,0,1)<birth.proportion){
			lucky.tip<-(sim.object$to[which(sim.object$tip)])[the.chosen.one]
			sim.object <- BirthSimObject(sim.object=sim.object, interval.length=interval.length, stop.time=depth.time, eps.constant=eps.param.indep, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc, max.tip.count=max.tip.count,lucky.tip=lucky.tip)
			max.tip.count<-max(sim.object$to)
			rate.track<-c(rate.track,turnover.splits[1])
			time.track<-c(time.track,depth.time)
			#Keeps track of the splits in the tree so we can add the last bit at the end if the starting root !== ending root:
			root.tracker<-rbind(root.tracker,c(sim.object$from[which(!sim.object$from%in%root.tracker[,1])][1],depth.time))
			phy <- sim2phylo(sim.object)
			if(verbose) {
				if(is.finite(max.time)){
					print(c(Ntip(phy),depth.time))
				}
			}
		}
		else {
			#If the lucky.tip is greater than the number of tips, then it is a death event and we rescale the value so that it actually corresponds to a tip:
			unlucky.tip <- (sim.object$to[which(sim.object$tip)])[the.chosen.one - alive]
			sim.object <- DeathSimObject(sim.object=sim.object, interval.length=interval.length, stop.time=depth.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc, unlucky.tip=unlucky.tip)
			if(!is.null(sim.object)){
				phy <- sim2phylo(sim.object) 
			}
			if (is.null(sim.object)) { 
				if(return.all.extinct) {
					return(NULL)
				}
				else {
					depth.time<-max.time
					root.tracker=matrix(c(1,depth.time),nrow=1,ncol=2)
					sim.object<-InitializeSimObject(birth, death, turnover.param.indep, eps.param.indep)
					phy <- sim2phylo(sim.object)
				}
			}
			if(verbose) {
				if(is.finite(max.time)){
					print(c(Ntip(phy),depth.time))
				}
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
	return(data.frame(from=c(1,1),to=c(2,3), edge.length=rep(0,2), tip=rep(TRUE,2), turnover.anc=birth, turnover.present=rep(birth-death,2), eps.anc=death, eps.present=rep(death/birth,2)))
}

GrowSimObject<-function(sim.object, interval.length) {
	sim.object$edge.length[which(sim.object$tip)]<-sim.object$edge.length[which(sim.object$tip)] + interval.length
	return(sim.object)
}

UpdateTerminalParams<-function(sim.object, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, split.times, turn.k, eps.k, turnover.splits, turnover.sigma.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent,  eps.splits, eps.sigma.anc) {
	for (descendant in sequence(dim(sim.object)[1])) {
		if(sim.object$tip[descendant]==TRUE){
			sim.object$turnover.present[descendant] <- SetParameter(stop.time, param.anc=sim.object$turnover.anc[descendant], sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits)
			sim.object$eps.present[descendant] <- SetParameter(stop.time, param.anc=sim.object$eps.anc[descendant], sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits)
		}
	}
	return(sim.object)
}

BirthSimObject<-function(sim.object, interval.length, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, split.times, turn.k, eps.k, turnover.splits, turnover.sigma.anc, eps.constant, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, eps.splits, eps.sigma.anc, max.tip.count, lucky.tip){
#	lucky.tip<-(sim.object$to[which(sim.object$tip)])[floor(runif(1,1,1+sum(sim.object$tip)))]
	lucky.edge<-which(sim.object$to==lucky.tip)
	sim.object<-GrowSimObject(sim.object, interval.length)
	sim.object<-UpdateTerminalParams(sim.object=sim.object, stop.time=stop.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
	sim.object<-rbind(sim.object, sim.object[lucky.edge,], sim.object[lucky.edge,])
	sim.object$tip[lucky.edge]<-FALSE
	#Update the ancestor for the new descendants:
	new.birth.anc <- GetNewAncParam(stop.time, parent.branch.length=sim.object$edge.length[lucky.edge], param.anc=sim.object$turnover.anc[lucky.edge], sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits, param.sigma.anc=turnover.sigma.anc)
	new.death.anc <- new.birth.anc * eps.constant
	descendants<-c(dim(sim.object)[1]-1, dim(sim.object)[1]) #the new taxa are in the last two rows, so just look at these
	for (descendant.index in sequence(length(descendants))) {
		descendant<-descendants[descendant.index]
		#sim.object$to[descendant]<-1+max(sim.object$to)
		sim.object$to[descendant]<-1+max.tip.count
		sim.object$from[descendant]<-lucky.tip
		sim.object$edge.length[descendant]<-0
		#Apply the ancestor of the new descendants:
		sim.object$turnover.anc[descendant]<-new.birth.anc
		sim.object$eps.anc[descendant]<-new.death.anc
		max.tip.count<-max(sim.object$to)
	}
	sim.object<-UpdateTerminalParams(sim.object=sim.object, stop.time=stop.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
	return(sim.object)
}

DeathSimObject<-function(sim.object, interval.length, stop.time, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, split.times, turn.k, eps.k, turnover.splits, turnover.sigma.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, eps.splits, eps.sigma.anc, unlucky.tip) {
	
	###########################REMAINING ISSUE###########################
	#Error is thrown when N<=2 during a death event -- not ideal. Although, 
	#we would have to live with this if we were still relying on ape...
	if(length(sim.object[,1])<=2){
		return(NULL)
	}
	#####################################################################
	else{
#		unlucky.tip<-(sim.object$to[which(sim.object$tip)])[floor(runif(1,1,1+sum(sim.object$tip)))]
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
			sim.object$to[which(sim.object$to==merged.edge)]<-surviving.descendant.of.merged.edge
			#Changes the new tip to a tip:
			sim.object$tip[which(sim.object$to==surviving.descendant.of.merged.edge)]<-TRUE
		}
		#If the sister is not a tip then the "to" must be changed for the ancestor of the non-tip sister
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
		sim.object<-UpdateTerminalParams(sim.object=sim.object, stop.time=stop.time, turnover.sigma.indep=turnover.sigma.indep, turnover.weight.anc=turnover.weight.anc, turnover.weight.logistic=turnover.weight.logistic, turnover.trend.exponent=turnover.trend.exponent, split.times=split.times, turn.k=turn.k, eps.k=eps.k, turnover.splits=turnover.splits, turnover.sigma.anc=turnover.sigma.anc, eps.sigma.indep=eps.sigma.indep, eps.weight.anc=eps.weight.anc, eps.weight.logistic=eps.weight.logistic, eps.trend.exponent=eps.trend.exponent, eps.splits=eps.splits, eps.sigma.anc=eps.sigma.anc)
		return(sim.object)
	}
}

GetNewAncParam<-function(stop.time, parent.branch.length, param.anc, sigma.indep, weight.anc, weight.logistic, trend.exponent, split.times, k, param.splits, param.sigma.anc){
	#A source of ambiguity comes from the odd rlnorm() function in R: it uses logmean and logsd. Here we assume that Rabosky was listing the values according to a normal and, therefore, we transformed these values to fit the desired values under a lognormal.
	s=sqrt(log((((param.sigma.anc*sqrt(parent.branch.length))/param.anc)^2)+1))
	m=log((param.anc^2)/sqrt(((param.sigma.anc*sqrt(parent.branch.length))^2)+(param.anc^2)))
	result<-rlnorm(1, m, s)
	#Here we attempt to mimic what Rabosky does. However, note that variance, if speciation rate is evolving in a Browian way, increases as sigma-squared*T, so sd should increase as sigma*sqrt(T). Based on the paper, Rabosky uses sigma*T. This is a problem: one expects that the variance at a given node based on the rate at the ancestral node should not vary based on how many intervening speciation events there are. Under the paramaterization of Rabosky, variance does change: the longer the branch, the more variance. So we do sd*sqrt(T)
	return(result)
}

SetParameter <- function(stop.time, param.anc, sigma.indep, weight.anc, weight.logistic, trend.exponent, split.times, k, param.splits) {
	position <- max(which(split.times>=stop.time),1)
	n.taxa <- 1+length(which(split.times>=stop.time))
	param.indep <- param.splits[position]
	param.starting <- (param.indep * (1 - weight.anc)) + (param.anc * weight.anc)
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
	birth.expected <- turnover / (1 - extinction.fraction)
	return(birth.expected)
}

SetDeath <- function(stop.time, turnover.param.anc, turnover.sigma.indep, turnover.weight.anc, turnover.weight.logistic, turnover.trend.exponent, eps.param.anc, eps.sigma.indep, eps.weight.anc, eps.weight.logistic, eps.trend.exponent, split.times, turn.k, eps.k, turnover.splits, eps.splits) {
	turnover<-SetParameter(stop.time=stop.time, param.anc=turnover.param.anc, sigma.indep=turnover.sigma.indep, weight.anc=turnover.weight.anc, weight.logistic=turnover.weight.logistic, trend.exponent=turnover.trend.exponent, split.times=split.times, k=turn.k, param.splits=turnover.splits)
	extinction.fraction<-SetParameter(stop.time=stop.time, param.anc=eps.param.anc, sigma.indep=eps.sigma.indep, weight.anc=eps.weight.anc, weight.logistic=eps.weight.logistic, trend.exponent=eps.trend.exponent, split.times=split.times, k=eps.k, param.splits=eps.splits)
	return(SetDeathGivenParameters(turnover, extinction.fraction))
}

SetDeathGivenParameters <- function(turnover, extinction.fraction) {
	death.expected<-(extinction.fraction * turnover) / (1 - extinction.fraction)
	return(death.expected)
}


