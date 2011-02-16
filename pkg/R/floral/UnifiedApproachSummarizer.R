library(diversitree) #obvious
library(sfsmisc) #for counting in binary
nchar=7

numberOfTransitionModelsForPartitionSize<-c(4,17,5+3*(2^3),5+3*(2^4),5+3*(2^5),5+3*(2^6),5+3*(2^7)) #make sure to include models that reduce to independent models: though they are present in smaller model sets, we might want independent transitions but dependent diversification
numberOfDiversificationModelsForPartitionSize<-c(19,12,6+3*(2^3),6+3*(2^4),6+3*(2^5),6+3*(2^6),6+3*(2^7))

summary.dataframe<-data.frame()


#three steps:
#Step 1
#	Get assignment of components: i.e., bisse for trait 1, musse for traits 1+2, musse for traits 3+5+6
#	For each component, try each of the possible models: i.e., for bisse, full, equal, uncorrelated; for 3 char models, maybe the "one is odd" models, split models, etc.
#   Record for each of these the model partition, combination type, model type (transition + diversification type), K, lnL, name of file containing the values
#
#Step 2
#	Get AIC of all the model combinations: i.e., AIC of aabbbcc where a=musse uniform, b=musse uncorrelated, etc.
#	Find the best model combinations
#	Report weights for different model types: how important is correlation between chars 2 and 3?
#
#Step 3
#	Get parameter estimates from best? equally good? average? models

#Step 1
#we will store these results in a dataframe
#char 1, 2, 3.... \t combosize\t model \t lnL \t K \t AIC \toutput file name
#but first have to make them
numberOfModels=0;
numberOfCompletedModels=0;
for (partitionIterator in 1:((2^nchar)-1)) { #do not want a model with no chars, so start with 1, not zero
	#print(digitsBase(partitionIterator,ndigits=nchar))
	partitionScheme<-digitsBase(partitionIterator,ndigits=nchar)
	partitionScheme<-partitionScheme[,1] #1 means include that char in partition, NA means do not
	partitionScheme<-partitionScheme*c(1:nchar) #so this will result in a vector with digits showing if the character is included
	partitionSize<-length(which(partitionScheme>0))
	partitionSchemeText=paste(partitionScheme,sep="",collapse="_")
	#mkdirCmd=paste("mkdir ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),sep="",collapse="")
	#suppressWarnings(system(mkdirCmd))
	for (transitionIndex in 1:numberOfTransitionModelsForPartitionSize[partitionSize]) { 
		for (diversificationIndex in 1:numberOfDiversificationModelsForPartitionSize[partitionSize]) {
			if (min(numberOfTransitionModelsForPartitionSize[partitionSize],numberOfDiversificationModelsForPartitionSize[partitionSize])>0) { #just to make sure we have models: for (1:1) and for (1:0) still run
				numberOfModels=numberOfModels+1
				#print(numberOfModels)
				nameRoot=paste("P",partitionSchemeText,"_T",transitionIndex,"_D",diversificationIndex,sep="",collapse="")
				#mkdirCmd=paste("mkdir ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,sep="",collapse="")
				#mkdirResult<-suppressWarnings(system(mkdirCmd))
				lsString=paste(paste("ls -1 ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,' | grep -c final.matrix.all',sep="",collapse=""))
				#print(lsString)
				finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
				#print(paste("finalMatrixAllCount = ",finalMatrixAllCount," for ",nameRoot))
				if(finalMatrixAllCount==1) {
					numberOfCompletedModels=numberOfCompletedModels+1
					lsString=paste(paste("ls -1 ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,' | grep final.matrix.all',sep="",collapse=""))
					finalMatrixName=suppressWarnings(system(lsString,intern=TRUE))
					suppressWarnings(rm(final.matrix.all)) #just to make sure anything we append is new
					suppressWarnings(rm(tmp.dataframe)) #ditto
					finalMatrixFullPath=paste(paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/',finalMatrixName,sep="")
					load(finalMatrixFullPath)
					dataFromFinalMatrix<-data.frame(final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),])
					names(dataFromFinalMatrix)<-c("lnLik","AIC","k_all","k_q","k_lambda","k_mu")
					tmp.dataframe<-data.frame(numberOfModels,numberOfCompletedModels,partitionSize,partitionSchemeText,transitionIndex,diversificationIndex)
					names(tmp.dataframe)<-c("modelID","numberCompletedModels","nChar","partition","T","D")
					tmp.dataframe<-cbind(tmp.dataframe,dataFromFinalMatrix)
					print(tmp.dataframe)
					summary.dataframe<-rbind(summary.dataframe,tmp.dataframe)
				}
			}
		}
	}
}
print(summary.dataframe)
write.table(summary.dataframe,file="Summary.txt",sep="\t")
save(summary.dataframe,file="Summary.Rsave",ascii=TRUE)

