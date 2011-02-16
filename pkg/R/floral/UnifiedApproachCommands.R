library(diversitree)
library(sfsmisc) 
library(ape)
library(partitions) #for converting from binary back to decimal

#the following 56 lines are truly stupid. However, otherwise I get an error thrown when i call make.cache.musse
source('../../../UnifiedApproachScripts/diversitree/R/asr-bisse.R')
source('../../../UnifiedApproachScripts/diversitree/R/asr-mkn.R')
source('../../../UnifiedApproachScripts/diversitree/R/asr-musse.R')
source('../../../UnifiedApproachScripts/diversitree/R/asr.R')
source('../../../UnifiedApproachScripts/diversitree/R/check.R')
source('../../../UnifiedApproachScripts/diversitree/R/clade-tree.R')
source('../../../UnifiedApproachScripts/diversitree/R/constrain.R')
source('../../../UnifiedApproachScripts/diversitree/R/diversitree-branches.R')
source('../../../UnifiedApproachScripts/diversitree/R/drop.tip.fixed.R')
source('../../../UnifiedApproachScripts/diversitree/R/history.R')
source('../../../UnifiedApproachScripts/diversitree/R/mcmc-norm.R')
source('../../../UnifiedApproachScripts/diversitree/R/mcmc-slice.R')
source('../../../UnifiedApproachScripts/diversitree/R/mcmc.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle-grid.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle-integer.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle-mixed.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle-optimize.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle-subplexR.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle-tgp.R')
source('../../../UnifiedApproachScripts/diversitree/R/mle.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bd-ode.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bd-split.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bd-t.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bd.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bisse-split.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bisse-t.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bisse-td.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bisse-unresolved.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bisse.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-bm.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-geosse.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-mkn.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-musse-split.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-musse-t.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-musse-td.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-musse.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-quasse-common.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-quasse-fftC.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-quasse-fftR.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-quasse-mol.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-quasse-split.R')
source('../../../UnifiedApproachScripts/diversitree/R/model-quasse.R')
source('../../../UnifiedApproachScripts/diversitree/R/plot-alt-extra.R')
source('../../../UnifiedApproachScripts/diversitree/R/plot-alt-util.R')
source('../../../UnifiedApproachScripts/diversitree/R/plot-alt.R')
source('../../../UnifiedApproachScripts/diversitree/R/profiles-plot.R')
source('../../../UnifiedApproachScripts/diversitree/R/simulate-bd.R')
source('../../../UnifiedApproachScripts/diversitree/R/simulate-bisse.R')
source('../../../UnifiedApproachScripts/diversitree/R/simulate-musse.R')
source('../../../UnifiedApproachScripts/diversitree/R/simulate-quasse.R')
source('../../../UnifiedApproachScripts/diversitree/R/simulation.R')
source('../../../UnifiedApproachScripts/diversitree/R/split-recycle.R')
source('../../../UnifiedApproachScripts/diversitree/R/split.R')
source('../../../UnifiedApproachScripts/diversitree/R/t.R')
source('../../../UnifiedApproachScripts/diversitree/R/td.R')
source('../../../UnifiedApproachScripts/diversitree/R/util.R')
source('../../../UnifiedApproachScripts/diversitree/R/')

print("Using command script v. Feb. 16 1:19 pm")

replaceextralist<-function(i) {
	#print("in replaceextralist")
	#print(paste("i = ",i))
	#print(paste("class i = ",class(i)," dim(i) =",dim(i), " length(i) = ",length(i)))
	assign("extralist",i,envir = .GlobalEnv)
}

make.musse.modifiedWithRootFixedAt1 <- function(tree, states, k, sampling.f=NULL, strict=FALSE,
                       safe=FALSE) { #we turn off strict so that we can still run even if we have just chars 1, 2, 3, 5, 6, 7, 8, for example
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  branches <- make.branches.musse(k, safe)
  root.p.vector=rep(0,k)
  root.p.vector[1]=1
  ll.musse <- function(pars, condition.surv=TRUE, root=ROOT.GIVEN,
                       root.p=root.p.vector, intermediates=FALSE) {
    if ( length(pars) != k*(k+1) )
      stop(sprintf("Invalid length parameters (expected %d)",
                   k*(k+1)))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    ll.xxsse(pars, cache, initial.conditions.musse, branches,
             condition.surv, root, root.p, intermediates)
  }

  ll <- function(pars, ...) ll.musse(pars, ...)
  class(ll) <- c("musse", "function")
  attr(ll, "k") <- k
  ll
}

doUnifiedRun<-function(P=P,T=T,D=D,S=partitionSize) {
	#first, make the data and tree files
	filename=prepData(P=P,T=T,D=D,S=S)

	#next, load the actual data
	data<-read.csv(file=paste(filename,".csv",sep=""),header=TRUE)
	states<-data[,2]
	names(states)<-data[,1]
	phy<-read.tree(file=paste(filename,"tree",".t",sep=""))
	
	#now start musse setup
	nAngiosperms=250000
	
	sampling.f<-rep(length(states)/nAngiosperms,2^S)
	results.vector.all<-c()
	lik <- make.musse.modifiedWithRootFixedAt1(tree=phy, states=states, k=2^S, sampling.f=sampling.f)
	#print("trying to assign extralist")
	assign("extralist",list(),envir = .GlobalEnv)
	#print(paste("first extralist is ",extralist))
	lik.trans <- modify_transitions(lik, type=T, S=S,extralist=extralist)
	#print(paste("second extralist is ",extralist))
	argnames(lik.trans)
	lik.final <- modify_diversification(lik.trans, type=D, S=S,extralist=extralist)
	#print(paste("third extralist is ",extralist))
	argnames(lik.final)
	p <- starting.point.musse(phy, 2^S)
	if (length(extralist)>0) {
		p <- starting.point.musse.extra(phy,2^S,argnames=argnames(lik.final)) #this is done as otherwise won't get right starting vector
	}
	fit.final <- find.mle(lik.final,p,method="subplex")
	save(fit.final, file=paste(filename,'.fit.final',sep=""), ascii=TRUE)
	print(fit.final)
	final.matrix<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),fit.final$par),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(fit.final$par))))
	save(final.matrix, file=paste(filename,'.final.matrix',sep=""), ascii=TRUE)
	rownames(final.matrix)<-paste("FINAL_",rownames(final.matrix),sep="") #to make it easier to grep
	print(formatC(final.matrix,format="f",digits=30,drop0trailing=TRUE))
	print(paste("FINAL_P ",P,sep=""))
	print(paste("FINAL_T ",T,sep=""))
	print(paste("FINAL_D ",D,sep=""))
	print(paste("FINAL_S ",S,sep=""))
	print(paste("FINAL_filename ",filename,sep=""))
	final.matrix.all<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),coef(fit.final,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(coef(fit.final,full=TRUE,extra=TRUE)))))
	save(final.matrix.all, file=paste(filename,'.final.matrix.all',sep=""), ascii=TRUE)
	rownames(final.matrix.all)<-paste("FINALALL_",rownames(final.matrix.all),sep="") #to make it easier to grep
	print(formatC(final.matrix.all,format="f",digits=30,drop0trailing=TRUE))
	
	
	#modify below this
	#diversificationVector<-getDiversificationRates(fit.final, type=diversificationType, partitionSize=S)
	#transitionVector<-getTransitionRates(fit.final,type=transitionType, partitionSize=S)
	#results.vector <- c(diversificationType, transitionType, logLik(fit.final), length(fit.final$par), getAIC(fit.final), diversificationVector, transitionVector)
	#save(results.vector,file=paste(filename,'raw.optim',sep=""),ascii=TRUE)
	#names(results.vector)<-c("diversificationType", "transitionType", "lnL", "K", "AIC", rep("div or trans ",length(results.vector)-5))
	#print(results.vector)
	#save(results.vector, file=paste(filename,'.optim',sep=""), ascii=TRUE)

}



modify_transitions<-function(lik=lik, type=1, S=S, extralist=extralist) {	
	print("in modify_transitions")
	print(paste("extralist is ",extralist))

	if (S==1) {
		if (type==1) {
			return(lik) #full_transition
		}
		else if (type==2) {
			return(constrain(lik, q12~q21)) #time-reversible
		}
		else if (type==3) {
			return(constrain(lik, q12~0)) #q01 is zero
		}
		else if (type==4) {
			return(constrain(lik, q21~0)) #q10 is zero
		}
		else {
			return("error")	
		}	
	}
	if (S==2) {
		if (type==1) {
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0)) #full_transition
		}
		else if (type==2) {
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q12~q34,q21~q43,q13~q24,q31~q42)) #uncorrelated_transition_rates
		}
		else if (type==3) {
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q12~q21, q13~q31, q24~q42, q34~q43)) #time-reversible
		}
		else if (type==4) {
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q12~q21, q13~q21, q24~q21, q31~q21, q34~q21, q42~q21, q43~q21)) #equal
		}
		else if (type==5) { #character 1 model two rate inflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q31~q21, q13~q12, q24~q12, q42~q12, q43~q12, q34~q12))
		}
		else if (type==6) { #character 1 model two rate outflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q13~q12, q31~q21, q24~q21, q42~q21, q43~q21, q34~q21))
		}
		else if (type==7) { #character 1 model two rate both unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q13~q12, q31~q12 ,q21~q12, q24~q42, q43~q42, q34~q42))
		}
		else if (type==8) { #character 2 model two rate inflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q12~q42, q21~q24, q43~q24, q34~q24, q31~q24, q13~q24))
		}
		else if (type==9) { #character 2 model two rate outflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q21~q24, q12~q42, q43~q42, q34~q42, q31~q42, q13~q42))
		}
		else if (type==10) { #character 2 model two rate both unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q21~q24, q12~q24 ,q42~q24, q43~q34, q31~q34, q13~q34))
		}
		else if (type==11) { #character 3 model two rate inflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q43~q13, q34~q31, q12~q31, q21~q31, q24~q31, q42~q31))
		}
		else if (type==12) { #character 3 model two rate outflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q34~q31, q43~q13, q12~q13, q21~q13, q24~q13, q42~q13))
		}
		else if (type==13) { #character 3 model two rate both unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q34~q31, q43~q31 ,q13~q31, q12~q21, q24~q21, q42~q21))
		}
		else if (type==14) { #character 4 model two rate inflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q24~q34, q42~q43, q31~q43, q13~q43, q12~q43, q21~q43))
		}
		else if (type==15) { #character 4 model two rate outflow unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q42~q43, q24~q34, q31~q34, q13~q34, q12~q34, q21~q34))
		}
		else if (type==16) { #character 4 model two rate both unique
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q42~q43, q24~q43 ,q34~q43, q31~q13, q12~q13, q21~q13))
		}
		else if (type==17) { #gain only 
			return(constrain(lik, q41~0,q14~0,q23~0,q32~0, q21~0, q31~0, q43~0, q42~0))
		}
		else {
			return("error")	
		}	
	}
	if (S>=3) {
		maxStringLength=nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
		#rather than typing manually all the restrictions, I will calculate this automatically
		constraintString="constrain(lik "
		for (charStateI in 1:((2^S))) { 
			for (charStateJ in (charStateI+1):((2^S))) { 
				if (charStateJ<=(2^S)) {
					binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
					binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
				#	print(binaryStateIVector)
				#	print(binaryStateJVector)
					numberMismatches=sum(1-(binaryStateIVector==binaryStateJVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
					if (numberMismatches>1) {
						constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~0, q',sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~0',sep="") 
					}
				}
			}
		}
		
		
		if (type==1) { #all transition rates the same
			for (charStateI in 1:((2^S))) { 
				for (charStateJ in (charStateI+1):((2^S))) { 
					if (charStateJ<=(2^S)) {
						binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
						binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
					#	print(binaryStateIVector)
					#	print(binaryStateJVector)
						numberMismatches=sum(1-(binaryStateIVector==binaryStateJVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
						if (numberMismatches==1) {
							constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~qAll, q',sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~qAll',sep="") 
						}
					}
				}
			}
			replaceextralist("qAll") #changes in parent frame
			constraintString=paste(constraintString,", extra='qAll'",sep="") 
		}
		else if (type==2) { #time reversible
			for (charStateI in 1:((2^S))) { 
				for (charStateJ in (charStateI+1):((2^S))) { 
					if (charStateJ<=(2^S)) {
						binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
						binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
					#	print(binaryStateIVector)
					#	print(binaryStateJVector)
						numberMismatches=sum(1-(binaryStateIVector==binaryStateJVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
						if (numberMismatches==1) {
							constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~q',sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sep="") 
						}
					}
				}
			}
			
		}
		else if (type==3) { #full, so nothing to do
		}
		else if (type>3 && type<=3+3*(2^S)) {
		
		
			modeltype=0; #1=inflow only different for one state, #2=outflow only different for one state, #3=outflow and inflow different for one state
			if (type<=3+(2^S)) { #do inflow model
				modeltype=1
				print("inflow rate for focal state unique")
			}
			else if (type<=3+2*(2^S)) { #do outflow model
				modeltype=2
				print("outflow rate for focal state unique")

			}
			else if (type<=3+3*(2^S)) { #do inflow and outflow model
				modeltype=3
				print("inflow and outflow rates for focal state unique")				
			}
			focalstate=1+(type-4)%%(2^S) #so if type is 4, it means use the first char, so (4-3)%%8=1. If type is 14, that means 10%%8 and so use the second car (but a different model, the outflow model)
		#	nonfocalratestateI=1
		#	nonfocalratestateJ=2
			#focalratestateI=1
			#focalratestateJ=focalstate
		#	if (focalstate==1) {
		#		nonfocalratestateI=3
		#	}
		#	else if (focalstate==2) {
		#		nonfocalratestateJ=3
		#	}
			for (charStateI in 1:((2^S))) { 
				for (charStateJ in 1:((2^S))) { #note we're starting from 1 here, too
					if (charStateJ<=(2^S)) {
						binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
						binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
						numberMismatches=sum(1-(binaryStateIVector==binaryStateJVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
						if (numberMismatches==1) {
					#		if (charStateI!=focalstate && charStateJ!=focalstate && charStateI!=nonfocalratestateI && charStateJ!=nonfocalratestateJ) { 
							if (charStateI!=focalstate && charStateJ!=focalstate) { 
								constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~qMost',sep="") 
							}
							else if (charStateI==focalstate) { #outflow rate
								if (modeltype!=1) { #so we have a unique outflow rate
									constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~qSpecial',sep="") 
								}
								else {
									constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~qMost',sep="") 
								}
							}
							else if (charStateJ==focalstate) { #inflow rate
								if (modeltype!=2) { #so we have a unique inflow rate
									constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~qSpecial',sep="") 
								}
								else {
									constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~qMost',sep="") 
								}
							}
						}
					}
				}
			}
			replaceextralist(c('qMost','qSpecial'))
			constraintString=paste(constraintString,", extra=c('qMost','qSpecial')",sep="") 
			
		}
		else if (type==(4+3*(2^S))) { #do gain-only model, otherwise full
			for (charStateI in 1:((2^S))) { 
				for (charStateJ in 1:((2^S))) { 
					if (charStateJ<=(2^S)) {
						binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
						binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
						numberMismatches=sum(1-(binaryStateIVector==binaryStateJVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
						if (numberMismatches==1) {
							if (sum(binaryStateJVector-binaryStateIVector)==-1) { #ie., 101 - 111 = -1, so it's been a loss
								constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~0',sep="") 
							}
						}
					}
				}
			}

		}
		else if (type==(5+3*(2^S))) { #do gain-only model, otherwise equal
			for (charStateI in 1:((2^S))) { 
				for (charStateJ in 1:((2^S))) { 
					if (charStateJ<=(2^S)) {
						binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
						binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
						numberMismatches=sum(1-(binaryStateIVector==binaryStateJVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
						if (numberMismatches==1) {
							if (sum(binaryStateJVector-binaryStateIVector)==-1) { #ie., 101 - 111 = -1, so it's been a loss
								constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~0',sep="") 
							}
							else {
								if (charStateI!=1 && charStateJ!=2) {
									constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~q',sprintf(paste("%0",maxStringLength,"d",sep=""),1),sprintf(paste("%0",maxStringLength,"d",sep=""),2),sep="") 
								}
							}
						}
					}
				}
			}
		}
		else {
			return("error")
		}
		constraintString=paste(constraintString,")",sep="") 
		print(paste("transition model: ",constraintString))
		return(eval(parse(text=constraintString)))
	}	
	 
}

modify_diversification<-function(lik=lik, type=1, S=S, extralist=extralist) {
	print(paste("diversification input extralist = ",extralist))
	maxStringLength=nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
	if (S==1) {
		if (type==1) {
			return(lik) #no change	
		}
		else if (type==2) {
			return(constrain(lik, lambda1 ~ lambda2)) #equal speciation
		}
		else if (type==3) {
			return(constrain(lik, mu1 ~ mu2)) #equal extinction
		}	
		else if (type==4) {
			return(constrain(lik, mu1~0, mu2~0)) #make_extinction_rates == 0
		}
		else if (type==5) {
			return(constrain(lik, lambda1 ~ lambda2, mu1 ~ mu2)) #one speciation rate, one extinction rate
		}
		else if (type==6) {
			return(constrain(lik, lambda1~lambda2, mu1~0)) #eqlam, mu1=0
		}
		else if (type==7) {
			return(constrain(lik, lambda1~lambda2, mu2~0)) #eqlam, mu2=0
		}
		else if (type==8){
			return(constrain(lik, lambda1~lambda2, mu1~0, mu2~0)) #eqlam, mu1=mu2=0
		}
		else if (type==9){
			return(constrain(lik, mu1~mu2, lambda1~0)) #eqmu, lam0=0
		}
		else if (type==10){
			return(constrain(lik, mu1~mu2, lambda2~0)) #eqmu, lam1=0
		}
		else if (type==11){
			return(constrain(lik, mu1~mu2, lambda1~0, lambda2~0)) #eqmu, lam0=lam1=0
		}
		else if (type==12){
			return(constrain(lik, lambda1~0, mu1~0)) #lam0=0, mu1=0, one way to get r0=0
		}
		else if (type==13){
			return(constrain(lik, lambda1~0, mu2~0)) #lam0=0, mu2=0, kind of silly
		}
		else if (type==14){
			return(constrain(lik, lambda1~0, mu1~0, mu2~0)) #only lam1 free
		}
		else if (type==15){
			return(constrain(lik, lambda2~0, mu1~0)) #lam1=0, mu1=0,kind of silly 
		}
		else if (type==16){
			return(constrain(lik, lambda2~0, mu2~0)) #lam1=0, mu2=0, one way to get r0=0
		}
		else if (type==17){
			return(constrain(lik, lambda2~0, mu1~0, mu2~0)) #only lam0 free
		}
		else if (type==18){
			return(constrain(lik, lambda1~0, lambda2~0, mu2~0)) #only mu1 free
		}
		else if (type==19){
			return(constrain(lik, lambda1~0, lambda2~0, mu2~0)) #only mu2 free
		}	
	}
	if (S==2) {
		if (type==1) {
			return(lik) #no change	
		}
		if (type==2) {
			return(constrain(lik, lambda1 ~ lambda2, lambda3 ~ lambda2, lambda4 ~ lambda2)) #equal speciation
		}
		if (type==3) {
			return(constrain(lik, mu1 ~ mu2, mu3 ~ mu2, mu4 ~ mu2)) #make_extinction_rates_equal
		}	
		if (type==4) {
			return(constrain(lik, mu1~0, mu2~0, mu3~0, mu4~0)) #make_extinction_rates = 0
		}
		if (type==5) {
			return(constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, lambda4 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1)) #one speciation rate, one extinction rate
		}
		if (type==6) {
			return(constrain(lik, lambda3 ~ lambda1, lambda2 ~ lambda1, lambda4 ~ lambda1,  mu1~0, mu2~0, mu3~0, mu4~0)) #one speciation rate, extinction rate=0
		}
		if (type==7) {
			return(constrain(lik, lambda1~lambda3, lambda2~lambda4, mu1~mu3, mu2~mu4)) #char1_not_affecting_diversification
		}
		if (type==8){
			return(constrain(lik, lambda1~lambda2, lambda3~lambda4, mu1~mu2, mu3~mu4)) #char2_not_affecting_diversification
		}
		if (type==9){
			return(constrain(lik, lambda2~lambda3, lambda4~lambda3, mu2~mu3, mu4~mu3)) #combo00_specially_affecting_diversification
		}
		if (type==10){
			return(constrain(lik, lambda1~lambda3, lambda4~lambda3, mu1~mu3, mu4~mu3)) #combo01_specially_affecting_diversification
		}
		if (type==11){
			return(constrain(lik, lambda1~lambda2, lambda4~lambda2, mu1~mu2, mu4~mu2)) #combo10_specially_affecting_diversification
		}
		if (type==12){
			return(constrain(lik, lambda1~lambda2, lambda3~lambda2, mu1~mu2, mu3~mu2)) #combo11_specially_affecting_diversification
		}
	}
	if (S>=3) {
		constraintString="constrain(lik "
		if (type==1) {
			#do nothing, is full model
		}
		else if (type==2) { #equal speciation
			for (charStateI in 2:((2^S))) { #we don't do char 1, as we set everything to that one
				constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~lambda',sprintf(paste("%0",maxStringLength,"d",sep=""),1),sep="") 
			}
		}
		else if (type==3) { #equal extinction
			for (charStateI in 2:((2^S))) { #we don't do char 1, as we set everything to that one
				constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~mu',sprintf(paste("%0",maxStringLength,"d",sep=""),1),sep="") 
			}
		}
		else if (type==4) { #one speciation rate, one extinction rate
			for (charStateI in 2:((2^S))) { #we don't do char 1, as we set everything to that one
				constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~mu',sprintf(paste("%0",maxStringLength,"d",sep=""),1),sep="") 
				constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~lambda',sprintf(paste("%0",maxStringLength,"d",sep=""),1),sep="") 
			}
		}
		else if (type==5) { #zero extinction, free speciation
			for (charStateI in 1:((2^S))) { 
				constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~0',sep="") 
			}
		}
		else if (type==6) { #zero extinction, equal speciation
			for (charStateI in 1:((2^S))) { 
				constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~0',sep="") 
			}
			for (charStateI in 2:((2^S))) { #we don't do char 1, as we set everything to that one
				constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~lambda',sprintf(paste("%0",maxStringLength,"d",sep=""),1),sep="") 
			}
		}
		else {
			modelIndex=type-7
			focalstate=1+modelIndex%%(2^S) #so it starts at 1
			modelFamily=1+floor(modelIndex/(2^S)) # 1, 2, 3....
			#modelFamily==1: one combo has unique speciation, all have equal extinction
			#modelFamily==2: one combo has unique speciation and extinction, all others have equal speciation rates and equal extinction rates
			#modelFamily==3: one combo has unique extinction, all have equal speciation
			nonfocalstate=2
			if (focalstate==2) {
				nonfocalstate=1
			}
			for (charStateI in 1:((2^S))) { 
				if (charStateI!=focalstate) {
					if (charStateI!=nonfocalstate) {
						constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~lambda',sprintf(paste("%0",maxStringLength,"d",sep=""),nonfocalstate),sep="") 
						constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~mu',sprintf(paste("%0",maxStringLength,"d",sep=""),nonfocalstate),sep="") 
					}
				}
				else {
					if (modelFamily==1) {
						constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~mu',sprintf(paste("%0",maxStringLength,"d",sep=""),nonfocalstate),sep="") 
					}
					if (modelFamily==3) {
						constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~lambda',sprintf(paste("%0",maxStringLength,"d",sep=""),nonfocalstate),sep="") 
					}
					#otherwise, these left free to vary
				}
			}
		}
		if (length(extralist)>0) {
			constraintString=paste(constraintString,", extra=c(",sep="") 		
			for (extraIndex in 1:length(extralist)) {
				constraintString=paste(constraintString,"'",extralist[extraIndex],"'",sep="")
				if (extraIndex<length(extralist)) {
					constraintString=paste(constraintString,", ",sep="")
				}
			}
			constraintString=paste(constraintString,")",sep="") 
		}
		constraintString=paste(constraintString,")",sep="") 
		print(paste("diversification model: ",constraintString))
		return(eval(parse(text=constraintString)))
	}
}


prepData<-function(P=P,T=T,D=D,S=S,sourcetraits="../../../SourceData/Steb7binaryJan19prunenoper_BCOPrune.csv") {
	partitionVector<-strsplit(P,split="_")
	print(partitionVector)
	partitionVector<-as.numeric(partitionVector[[1]])
	print(partitionVector)
	charsToInclude<-partitionVector[which(partitionVector>0)]
	stopifnot(S==length(charsToInclude))
	file<-sourcetraits
	phy<-"../../../SourceData/floral_1.nex"
	tree<-read.nexus(phy)
	data<-read.csv(file)
	colnamesVector<-colnames(data)
	names<-colnamesVector[2:length(colnamesVector)]
	print(names)
	#keep only 1 (1st 25%) and 2 (2nd 25%) in Selection column
	#subdata<-subset(data,data$Selection<=2)
	#already done, so just
	subdata<-data #to minimize recoding
	
	#delete taxa with missing data anywhere
	for (colToExamine in 2:length(colnamesVector)) {
		subdata<-subdata[subdata[,colToExamine]!="?",]
	}
	
	#do recoding. Put chars at end of vector
	subdata[,(dim(subdata)[2]+1)]=subdata[,(charsToInclude[1]+1) ]
	if (length(charsToInclude)>=2) {
		for (charIndex in 2:length(charsToInclude)) {
			subdata[,(dim(subdata)[2])]=paste(subdata[,(dim(subdata)[2])], subdata[,(charsToInclude[charIndex]+1) ],sep="")
		}
	}
	#now make a new column that takes the 0, or 0 1, or 0 1 0 1, etc. and maps them into 1, 2, 3, 4
	#subdata[,(dim(subdata)[2]+1)]=1+todec(as.vector(unlist(strsplit(as.character(unlist(subdata[,(dim(subdata)[2]) ])," "))),mode="numeric"))
	subdata[,(dim(subdata)[2]+1)]<-NA
	#print(dim(subdata)[1])
	for (rowIndex in 1:dim(subdata)[1]) {
		rawData<-subdata[rowIndex,(dim(subdata)[2]-1)]
		rawData<-as.character(rawData)
		rawDataSplit<-strsplit(rawData,"")[[1]]
		rawDataSplit<-as.numeric(rawDataSplit)
		subdata[rowIndex,(dim(subdata)[2])]<-1+todec(rawDataSplit)
	}
	#if there are NA for the new char, this will take them out
	#if not, it will do nothing
	presentTaxa<-as.vector(subdata$Name_in_tree)
	to.drop <- setdiff(tree$tip.label, presentTaxa)
	tree2<-drop.tip(tree,to.drop)
	char<-subdata[,(dim(subdata)[2])]
	names(char)<-subdata$Name_in_tree
	colnamesVector<-colnames(data)
	print(colnamesVector)
	finalname=colnamesVector[(charsToInclude[1]+1)]
	print(finalname)
	if (length(charsToInclude)>=2) {
		for (charIndex in 2:length(charsToInclude)) {
			finalname=paste(finalname, colnamesVector[(charsToInclude[charIndex]+1) ],sep="_")
		}
	}
	#write out datafile for the character with matching tree
	write.csv(char,file=paste(finalname,".csv",sep=""))
	write.tree(tree2,file=paste(finalname,"tree",".t",sep=""))
	return(finalname)
}

