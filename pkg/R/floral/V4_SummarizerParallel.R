library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
source("V4_UtilityFns.R")
library(foreach)


focalVectorList<-getAllInterestingFocalVectorsStringsEfficient(S)

tVector=c(1:dim(transitionModels)[2])
dVector=c(1:dim(diversificationModels)[2])


summarizeIndiv<-function(actualT,actualD,focalVectorList) {
	
	loadedOld<-FALSE
	totalRuns<-0
	completedRuns<-0
	maxStringLength=nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
	
	
	summary.dataframe<-data.frame()
	
	suppressWarnings(system("mkdir -p ../Summaries"))
	for (focalIndex in 1:length(focalVectorList)) {
		focalVector<-stringToVector(unlist(focalVectorList[[focalIndex]]))
		for (transitionModelIndex in actualT:actualT) {
			for (diversificationModelIndex in actualD:actualD) {
				if (numberFocalCombos(focalVector) >= transitionModels$min_focalcombos[transitionModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
					if(numberFocalCombos(focalVector) >= diversificationModels$min_focalcombos[diversificationModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
						totalRuns<-totalRuns+1
						#yay! Now we can run!
						tryLoad<-TRUE
						if (loadedOld==TRUE) { #see if we've already loaded this
							if(length(which(old.summary.dataframe$T==transitionModelIndex & old.summary.dataframe$D==diversificationModelIndex & old.summary.dataframe$focal==paste(getFocalSummaryLabel(focalVector,S=7,any="x"),sep="",collapse="") ))==1) {
								tryLoad<-FALSE #it's already in the old.summary.dataframe
								completedRuns<-completedRuns+1
								print(paste("already have completed run ",completedRuns,"/",totalRuns,sep=""))
							}
						}
						if(tryLoad==TRUE) {
							nameRoot<-paste("T",transitionModelIndex,"_D",diversificationModelIndex,"_",vectorToString(getFocalSummaryLabel(focalVector,S,"x")),sep="",collapse="")
							dirRoot<-paste("../ActualRuns/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,"/",nameRoot,sep="",collapse="")
							lsString=paste(paste("ls -1 ",dirRoot,' | grep -c final.matrix.all',sep="",collapse=""))
							print(lsString)
							finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
							if(finalMatrixAllCount>0) {
								completedRuns<-completedRuns+1
								suppressWarnings(rm(final.matrix.all)) #just to make sure anything we append is new
								suppressWarnings(rm(tmp.dataframe)) #ditto
								finalMatrixFullPath<-paste(dirRoot,"/Steb1Perianth_Steb2PerFusSDS_Steb3SymSDS_Steb4StamNo_Steb5Syncarpy_Steb6SeedNo_Steb8Ovary.final.matrix.all",sep="")
								load(finalMatrixFullPath)
								qIndices<-grep("^q\\d",row.names(final.matrix.all),perl=TRUE)
								lambdaIndices<-grep("^lambda\\d",row.names(final.matrix.all),perl=TRUE)
								muIndices<-grep("^mu\\d",row.names(final.matrix.all),perl=TRUE)
								tmp.dataframe<-data.frame(paste(getFocalSummaryLabel(focalVector,S=7,any="x"),sep="",collapse=""),transitionModelIndex,transitionModels[transitionModelIndex,4],diversificationModelIndex,diversificationModels[diversificationModelIndex,5],final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),1],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),1]) 
								names(tmp.dataframe)<-c("focal","T","TransitionModel","D","DiversificationModel","lnLik","AIC","k_all","k_q","k_lambda","k_mu")	
								tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(final.matrix.all[qIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[qIndices,1])))),data.frame(matrix(final.matrix.all[lambdaIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[lambdaIndices,1])))),data.frame(matrix(final.matrix.all[muIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[muIndices,1])))))
								summary.dataframe<-rbind(summary.dataframe,tmp.dataframe)
								print(paste("loaded completed run ",completedRuns,"/",totalRuns,sep=""))
							}
						}
					}
				}
			}
		}
	}
	
	summary.dataframe$focal<-as.character(summary.dataframe$focal) #because we do not want a factor
	deltaAIC<-summary.dataframe$AIC-min(summary.dataframe$AIC)
	relativeLikelihood<-exp(-0.5 * deltaAIC)
	AICweight<-relativeLikelihood/sum(relativeLikelihood)
	summary.dataframe<-cbind(deltaAIC,AICweight,summary.dataframe)
	
	#now time to make nice names for things
	for (charStateI in 1:((2^S))) { 
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		iLabelShort<-sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI)
		iLabelLong<-vectorToString(binaryStateIVector)
		names(summary.dataframe)[which(names(summary.dataframe) == paste("lambda",iLabelShort,sep="",collapse=""))]<-paste("lambda",iLabelLong,sep="",collapse="")
		names(summary.dataframe)[which(names(summary.dataframe) == paste("mu",iLabelShort,sep="",collapse=""))]<-paste("mu",iLabelLong,sep="",collapse="")
		print(paste("changing names for ",iLabelLong))
		for (charStateJ in 1:((2^S))) { 
			binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
			jLabelShort<-sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ)
			jLabelLong<-vectorToString(binaryStateJVector)
			numberMismatches=vectorMismatch(binaryStateIVector,binaryStateJVector) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
			if (numberMismatches==1) {
				names(summary.dataframe)[which(names(summary.dataframe) == paste("q",iLabelShort,jLabelShort,sep="",collapse=""))]<-paste("q",iLabelLong,"_",jLabelLong,sep="",collapse="")			
			}
			else {
				names(summary.dataframe)[which(names(summary.dataframe) == paste("q",iLabelShort,jLabelShort,sep="",collapse=""))]<-paste("q",iLabelLong,"_",jLabelLong,"_disallowed",sep="",collapse="")			
			}
		}
	}
	
	save(summary.dataframe,file=paste("../Summaries/RateSummaryT",actualT,"D",actualD,".Rsave",compress=TRUE)
	return(paste("T",actualT,"D",actualD,completedRuns,"/",totalRuns))
}

finalresult<-foreach(actualT=tVector) %:% foreach(actualD=dVector) %doPar% { summarizeIndiv(actualT,actualD,focalVectorList) }
print(finalResult)
