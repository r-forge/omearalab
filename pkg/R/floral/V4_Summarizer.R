library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
source("V4_UtilityFns.R")


focalVectorList<-getAllInterestingFocalVectorsStringsEfficient(S)

totalRuns<-0
completedRuns<-0

summary.dataframe<-data.frame()

suppressWarnings(system("mkdir -p ../Summaries"))
for (focalIndex in 1:length(focalVectorList)) {
	focalVector<-stringToVector(unlist(focalVectorList[[focalIndex]]))
	for (transitionModelIndex in 1:dim(transitionModels)[1]) {
		for (diversificationModelIndex in 1:dim(diversificationModels)[1]) {
			if (numberFocalCombos(focalVector) >= transitionModels$min_focalcombos[transitionModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
				if(numberFocalCombos(focalVector) >= diversificationModels$min_focalcombos[diversificationModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
					totalRuns<-totalRuns+1
					#yay! Now we can run!
					nameRoot<-paste("T",transitionModelIndex,"_D",diversificationModelIndex,"_",vectorToString(getFocalSummaryLabel(focalVector,S,"x")),sep="",collapse="")
					dirRoot<-paste("../ActualRuns/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,"/",nameRoot,sep="",collapse="")
					lsString=paste(paste("ls -1 ",dirRoot,' | grep -c final.matrix.all',sep="",collapse=""))
					print(lsString)
					finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
					maxStringLength=nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
					if(finalMatrixAllCount>0) {
						completedRuns<-completedRuns+1
						suppressWarnings(rm(final.matrix.all)) #just to make sure anything we append is new
						suppressWarnings(rm(tmp.dataframe)) #ditto
						finalMatrixFullPath<-paste(dirRoot,"/Steb1Perianth_Steb2PerFusSDS_Steb3SymSDS_Steb4StamNo_Steb5Syncarpy_Steb6SeedNo_Steb8Ovary.final.matrix.all",sep="")
						load(finalMatrixFullPath)
						qIndices<-grep("^q\\d",row.names(final.matrix.all),perl=TRUE)
						lambdaIndices<-grep("^lambda\\d",row.names(final.matrix.all),perl=TRUE)
						muIndices<-grep("^mu\\d",row.names(final.matrix.all),perl=TRUE)
						tmp.dataframe<-data.frame(vectorToString(focalVector),transitionModelIndex,transitionModels[transitionModelIndex,4],diversificationModelIndex,diversificationModels[diversificationModelIndex,5],final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),1],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),1]) 
						names(tmp.dataframe)<-c("focal","T","TransitionModel","D","DiversificationModel","lnLik","AIC","k_all","k_q","k_lambda","k_mu")	
						tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(final.matrix.all[qIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[qIndices,1])))),data.frame(matrix(final.matrix.all[lambdaIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[lambdaIndices,1])))),data.frame(matrix(final.matrix.all[muIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[muIndices,1])))))
						summary.dataframe<-rbind(summary.dataframe,tmp.dataframe)
						print(paste("completed run ",completedRuns,"/",totalRuns,sep=""))
						write.table(summary.dataframe,file="../Summaries/RateSummary.txt",sep="\t")
						save(summary.dataframe,file="../Summaries/RateSummary.Rsave",compress=TRUE)
					}
				}
			}
		}
	}
}
