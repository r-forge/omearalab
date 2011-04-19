rm(list=ls())
library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions) #for converting from binary back to decimal
nchar=7

numberOfTransitionModelsForPartitionSize<-c(4,17,5+3*(2^3),5+3*(2^4),5+3*(2^5),5+3*(2^6),5+3*(2^7)) #make sure to include models that reduce to independent models: though they are present in smaller model sets, we might want independent transitions but dependent diversification
numberOfDiversificationModelsForPartitionSize<-c(19,12,6+3*(2^3),6+3*(2^4),6+3*(2^5),6+3*(2^6),6+3*(2^7))

summary.dataframe<-data.frame()

#we're going to convert the final rate estimates (i.e., rate of going from 01***** to 11***** for model 1_2_0_0_0_0_0) 
#   into the rates on the full 128 x 128 matrix: rate 01***** -> 11***** is used for rate of 0111111 -> 1111111, and 
#   for 0100001 -> 1100001, and so forth. The problem is this has no information on the rate of 0000000 -> 0000001, as
#   these changes aren't seen. But that may make the rates more combinable across partitions, always a good thing
#we just store raw AIC, not ∆AIC, as that will be based on summing AIC from models with different characters

numberOfModels=0;
numberOfCompletedModels=0;
for (partitionIterator in ((2^nchar)-1):((2^nchar)-1)) { #do not want a model with no chars, so start with 1, not zero
	#print(digitsBase(partitionIterator,ndigits=nchar))
	partitionScheme<-digitsBase(partitionIterator,ndigits=nchar)
	partitionScheme<-partitionScheme[,1] #1 means include that char in partition, NA means do not
	partitionScheme<-partitionScheme*c(1:nchar) #so this will result in a vector with digits showing if the character is included
	partitionSize<-length(which(partitionScheme>0))
	partitionSchemeText=paste(partitionScheme,sep="",collapse="_")
	print(partitionSchemeText)
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
				lsString=paste(paste("ls -1 ",paste("/Users/bomeara/Sites/RunsFebruary2011/ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,' | grep -c final.matrix.all',sep="",collapse=""))
				#print(lsString)
				finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
				#print(paste("finalMatrixAllCount = ",finalMatrixAllCount," for ",nameRoot))
				if(finalMatrixAllCount==1) {
					rateMatrix<-matrix(nrow=2^nchar, ncol=2^nchar) #we'll later make this a vector and stick it with the other info as a row in the data.frame. Of course most of these entries are forced to be zero, so we could be more efficient in storage by dropping them. But my efficiency is worth more
					#fill in the obligately zero rates
					for (charStateRow in 1:(2^nchar)) { 
						for (charStateCol in 1:(2^nchar)) { 
							binaryStateRowVector<-digitsBase(charStateRow-1,ndigits=nchar)[,1]
							binaryStateColVector<-digitsBase(charStateCol-1,ndigits=nchar)[,1]
							numberMismatches=sum(1-(binaryStateRowVector==binaryStateColVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
							if (numberMismatches>1) {
								rateMatrix[charStateRow,charStateCol]=0 #cause no leaps across multiple steps
							}
							else if (numberMismatches==0) {
								rateMatrix[charStateRow,charStateCol]=-1 #it's obviously not -1, but this will let us flag the diagonals later
							}
							#else it's a character combo differing by one change
						}
					}
					birthVector<-rep(NA,2^nchar)
					deathVector<-rep(NA,2^nchar)
					numberOfCompletedModels=numberOfCompletedModels+1
					lsString=paste(paste("ls -1 ",paste("/Users/bomeara/Sites/RunsFebruary2011/ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,' | grep final.matrix.all',sep="",collapse=""))
					finalMatrixName=suppressWarnings(system(lsString,intern=TRUE))
					suppressWarnings(rm(final.matrix.all)) #just to make sure anything we append is new
					suppressWarnings(rm(tmp.dataframe)) #ditto
					finalMatrixFullPath=paste(paste("/Users/bomeara/Sites/RunsFebruary2011/ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/',finalMatrixName,sep="")
					load(finalMatrixFullPath)
					#dataFromFinalMatrix<-data.frame(final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),1],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),1])
					#names(dataFromFinalMatrix)<-c("lnLik","AIC","k_all","k_q","k_lambda","k_mu")
					#Now the fun bits: getting the rate estimates
					qIndices<-grep("^q",row.names(final.matrix.all))
					lambdaIndices<-grep("^lambda",row.names(final.matrix.all))
					muIndices<-grep("^mu",row.names(final.matrix.all))
					#have to convert stuff like q13 into q0*0****->1*0****
					maxStringLength=nchar(2^partitionSize) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.

					for (charStateRow in 1:(2^nchar)) { 
						binaryStateRowVector<-digitsBase(charStateRow-1,ndigits=nchar)[,1]
						relevantNumbers<-binaryStateRowVector[which(partitionScheme>0)]
						rowStateLabel<-sprintf(paste("%0",maxStringLength,"d",sep=""),1+todec(relevantNumbers))
						birthVector[charStateRow]=final.matrix.all[which(row.names(final.matrix.all)==paste("lambda",rowStateLabel,sep="")),1]
						names(birthVector)[charStateRow]=paste("lambda",paste(binaryStateRowVector,sep="",collapse=""),sep="",collapse="")
						deathVector[charStateRow]=final.matrix.all[which(row.names(final.matrix.all)==paste("mu",rowStateLabel,sep="")),1]
						names(deathVector)[charStateRow]=paste("mu",paste(binaryStateRowVector,sep="",collapse=""),sep="",collapse="")
						for (charStateCol in 1:(2^nchar)) { 
							binaryStateColVector<-digitsBase(charStateCol-1,ndigits=nchar)[,1]
							relevantNumbers<-binaryStateColVector[which(partitionScheme>0)]
							colStateLabel<-sprintf(paste("%0",maxStringLength,"d",sep=""),1+todec(relevantNumbers))					
							numberMismatchesTotal=sum(1-(binaryStateRowVector==binaryStateColVector)) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
							numberMismatchesRelevantSites=sum(1-(binaryStateRowVector[which(partitionScheme>0)]==binaryStateColVector[which(partitionScheme>0)]))
							if (numberMismatchesTotal==1) {
								if (numberMismatchesRelevantSites==1 ) {
									rateMatrix[charStateRow,charStateCol] = final.matrix.all[which(row.names(final.matrix.all)==paste("q",rowStateLabel,colStateLabel,sep="")),1]#row = from,col=to
								}
							}
						}
					}	
					tmp.dataframe<-data.frame(final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),1],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),1],numberOfModels,numberOfCompletedModels,partitionSize,partitionSchemeText,transitionIndex,diversificationIndex) 
					names(tmp.dataframe)<-c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","modelID","numberCompletedModels","nChar","partition","T","D")	
					tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(birthVector,nrow=1,dimnames=list("",names(birthVector)))))
					tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(deathVector,nrow=1,dimnames=list("",names(deathVector)))))
					for (charStateRow in 1:(2^nchar)) { 
						binaryStateRowVector<-digitsBase(charStateRow-1,ndigits=nchar)[,1]
						ratevectorThisRow<-rateMatrix[charStateRow,]
						ratevectorThisRow[charStateRow]<-0-sum(ratevectorThisRow[-1*charStateRow],na.rm=TRUE) #so the diagonal is the negative of the sum of the others. 
						for (charStateCol in 1:(2^nchar)) { 
							binaryStateColVector<-digitsBase(charStateCol-1,ndigits=nchar)[,1]
							names(ratevectorThisRow)[charStateCol]<-paste("q",paste(binaryStateRowVector,sep="",collapse=""),"_",paste(binaryStateColVector,sep="",collapse=""),sep="",collapse="")
						}
						tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(ratevectorThisRow,nrow=1,dimnames=list("",names(ratevectorThisRow)))))
					}
					#print(tmp.dataframe)
					summary.dataframe<-rbind(summary.dataframe,tmp.dataframe)
					#write.table(summary.dataframe,file="RateSummary.txt",sep="\t")
					#save(summary.dataframe,file="RateSummary.Rsave",ascii=TRUE)

				}
			}
		}
	}
	#write.table(summary.dataframe,file="../Summaries/RateSummary.txt",sep="\t")
	#save(summary.dataframe,file="../Summaries/RateSummary.Rsave",ascii=TRUE)

}

print(summary.dataframe)
write.table(summary.dataframe,file="../Summaries/RateSummary.txt",sep="\t")
save(summary.dataframe,file="../Summaries/RateSummary.Rsave",compress=TRUE)

