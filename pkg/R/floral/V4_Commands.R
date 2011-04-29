library(diversitree)
library(sfsmisc) 
library(ape)
library(partitions) #for converting from binary back to decimal
source("V4_UtilityFns.R")



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

print("Using command script v. April 28 1:48 pm")

replaceextralist<-function(i) {
	#print("in replaceextralist")
	#print(paste("i = ",i))
	#print(paste("class i = ",class(i)," dim(i) =",dim(i), " length(i) = ",length(i)))
	assign("extralist",i,envir = .GlobalEnv)
}

make.musse.modifiedWithRootFixedAt1 <- function(tree, states, k, sampling.f=NULL, strict=FALSE,
                       safe=FALSE) { #we turn off strict so that we can still run even if we have just chars 1, 2, 3, 5, 6, 7, 8, for example. And remember that states start with 1, not 0, for musse.
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

#F=focal states
doUnifiedRun<-function(F=F, T=T,D=D,S=partitionSize) {
	#first, make the data and tree files
	filename=prepData(P="1_2_3_4_5_6_7",F=F,T=T,D=D,S=S)
	focalVector=F
	if (length(focalVector==1)) {
		focalVector<-stringToVector(focalVector)
	}
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
	assign("extralist",list(),envir = .GlobalEnv) #naughty
	#print(paste("first extralist is ",extralist))
	lik.trans <- modify_transitions(lik, type=T, F=F, S=S,extralist=extralist)
	#print(paste("second extralist is ",extralist))
	argnames(lik.trans)
	lik.final <- modify_diversification(lik.trans, type=D, F=F, S=S,extralist=extralist)
	#print(paste("third extralist is ",extralist))
	argnames(lik.final)
	p <- starting.point.musse(phy, 2^S)
	if (length(extralist)>0) {
		p <- starting.point.musse.extra(phy,2^S,argnames=argnames(lik.final)) #this is done as otherwise won't get right starting vector
	}
	fit.final <- find.mle(lik.final,p,method="subplex")
	#save(fit.final, file=paste(filename,'.fit.final',sep=""), compress=TRUE)
	print(fit.final)
	final.matrix<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),fit.final$par),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(fit.final$par))))
	#save(final.matrix, file=paste(filename,'.final.matrix',sep=""), compress=TRUE)
	rownames(final.matrix)<-paste("FINAL_",rownames(final.matrix),sep="") #to make it easier to grep
	print(formatC(final.matrix,format="f",digits=30,drop0trailing=TRUE))
	print(paste("FINAL_F ",F,sep=""))
	print(paste("FINAL_T ",T,sep=""))
	print(paste("FINAL_D ",D,sep=""))
	print(paste("FINAL_S ",S,sep=""))
	print(paste("FINAL_filename ",filename,sep=""))
	final.matrix.all<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),coef(fit.final,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(coef(fit.final,full=TRUE,extra=TRUE)))))
	save(final.matrix.all, file=paste(filename,'.final.matrix.all',sep=""), compress=TRUE)
	#rownames(final.matrix.all)<-paste("FINALALL_",rownames(final.matrix.all),sep="") #to make it easier to grep
	#print(formatC(final.matrix.all,format="f",digits=30,drop0trailing=TRUE))
	
	
	#modify below this
	#diversificationVector<-getDiversificationRates(fit.final, type=diversificationType, partitionSize=S)
	#transitionVector<-getTransitionRates(fit.final,type=transitionType, partitionSize=S)
	#results.vector <- c(diversificationType, transitionType, logLik(fit.final), length(fit.final$par), getAIC(fit.final), diversificationVector, transitionVector)
	#save(results.vector,file=paste(filename,'raw.optim',sep=""),ascii=TRUE)
	#names(results.vector)<-c("diversificationType", "transitionType", "lnL", "K", "AIC", rep("div or trans ",length(results.vector)-5))
	#print(results.vector)
	#save(results.vector, file=paste(filename,'.optim',sep=""), ascii=TRUE)

}


modify_transitions<-function(lik=lik, type=1, F=F, S=S, extralist=extralist) {
	#remember to see transition models in V*_UtilityFns.R
	maxStringLength=nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.



	#rather than typing manually all the restrictions, I will calculate this automatically
	constraintString="constrain(lik "
	for (charStateI in 1:((2^S))) { 
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		for (charStateJ in (charStateI+1):((2^S))) { 
			if (charStateJ<=(2^S)) {
				binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
				numberMismatches=vectorMismatch(binaryStateIVector,binaryStateJVector) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
				if (numberMismatches>1) {
					constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~0, q',sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~0',sep="") 
				}
			}
		}
	}
	
	#now for the actual model
	for (charStateI in 1:((2^S))) { 
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		fromFocal<-FALSE
		if( focalVector[comboAsDecimal(binaryStateIVector,S)] == 1) {
			fromFocal<-TRUE
		}
		for (charStateJ in 1:((2^S))) { #note we're starting from 1 here, too
			binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
			numberMismatches=vectorMismatch(binaryStateIVector,binaryStateJVector) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
			if (numberMismatches==1) {
				toFocal<-FALSE
				if( focalVector[comboAsDecimal(binaryStateJVector,S)] == 1) {
					toFocal<-TRUE
				}
				qString<-"~qMost"
				if(fromFocal) {
					if (toFocal) {
						qString<-qFFbyModel[type] #see utility file
					}
					else {
						qString<-qFNbyModel[type]
					}
				}
				else {
					if (toFocal) {
						qString<-qNFbyModel[type]
					}
					else {
						qString<-qNNbyModel[type]
					}
				}
				constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),qString,sep="") 
			}
		}
	}
	constraintString=paste(constraintString,")",sep="") 
	print(paste("transition model: ",constraintString))
	return(eval(parse(text=constraintString)))
}

modify_diversification<-function(lik=lik, type=1, F=F, S=S, extralist=extralist) {
	print(paste("diversification input extralist = ",extralist))
	maxStringLength<-nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
	constraintString<-"constrain(lik "
	for (charStateI in 1:((2^S))) { 
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		isFocal<-FALSE
		if( focalVector[comboAsDecimal(binaryStateIVector,S)] == 1) {
			isFocal<-TRUE
		}
		if (isFocal) {
			constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),bFbyModel[type],sep="") 
			constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),dFbyModel[type],sep="") 
		}
		else {
			constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),bNbyModel[type],sep="") 
			constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),dNbyModel[type],sep="") 
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

prepData<-function(P=P,F=F,T=T,D=D,S=S,sourcetraits="../../../SourceData/Steb7binaryJan19prunenoper_BCOPrune.csv") {
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

