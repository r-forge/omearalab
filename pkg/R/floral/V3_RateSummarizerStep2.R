rm(list=ls())
library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions) #for converting from binary back to decimal
nchar=7

load(file="../Summaries/RateSummary.Rsave") #welcome, summary.dataframe
numberOfTransitionModelsForPartitionSize<-c(4,17,5+3*(2^3),5+3*(2^4),5+3*(2^5),5+3*(2^6),5+3*(2^7)) #make sure to include models that reduce to independent models: though they are present in smaller model sets, we might want independent transitions but dependent diversification
numberOfDiversificationModelsForPartitionSize<-c(19,12,6+3*(2^3),6+3*(2^4),6+3*(2^5),6+3*(2^6),6+3*(2^7))

summary2.dataframe<-data.frame()

qNamesVector=c()
for (charStateRow in 1:(2^nchar)) { 
	for (charStateCol in 1:(2^nchar)) { 
		binaryStateRowVector<-digitsBase(charStateRow-1,ndigits=nchar)[,1]
		binaryStateColVector<-digitsBase(charStateCol-1,ndigits=nchar)[,1]
		numberMismatches=sum(1-(binaryStateRowVector==binaryStateColVector))
		if (numberMismatches<=1) { #returns diagonal and possible transitions
			qNamesVector<-append(qNamesVector,paste("q",paste(binaryStateRowVector,sep="",collapse=""),"_",paste(binaryStateColVector,sep="",collapse=""),sep="",collapse=""))
		}
	}
}

for (partitionIterator in ((2^nchar)-1):((2^nchar)-1)) { #do not want a model with no chars, so start with 1, not zero
	#print(digitsBase(partitionIterator,ndigits=nchar))
	partitionScheme<-digitsBase(partitionIterator,ndigits=nchar)
	partitionScheme<-partitionScheme[,1] #1 means include that char in partition, NA means do not
	partitionScheme<-partitionScheme*c(1:nchar) #so this will result in a vector with digits showing if the character is included
	partitionSize<-length(which(partitionScheme>0))
	partitionSchemeText=paste(partitionScheme,sep="",collapse="_")
	print(partitionSchemeText)
	subset.summary.dataframe<-subset(summary.dataframe,summary.dataframe$partition==partitionSchemeText)
	if (dim(subset.summary.dataframe)[1]>0) {
		deltaAIC=subset.summary.dataframe$AIC-min(subset.summary.dataframe$AIC,na.rm=TRUE)
		relativeLikelihood=exp(-0.5 * deltaAIC)
		AICweight=relativeLikelihood/sum(relativeLikelihood,na.rm=TRUE)
		subset.summary.dataframe<-cbind(subset.summary.dataframe,deltaAIC,AICweight)
		write.table(subset.summary.dataframe,file=paste("../Summaries/IndivRateSummary_",partitionSchemeText,".txt",sep=""),sep="\t")
		#save(subset.summary.dataframe,file=paste("../Summaries/IndivRateSummary_",partitionSchemeText,".Rsave",sep=""),ascii=TRUE)
		allBirthVectors<-subset.summary.dataframe[,grep("^lambda\\d",names(subset.summary.dataframe))]
		allDeathVectors<-subset.summary.dataframe[,grep("^mu\\d",names(subset.summary.dataframe))]
		allTransitionVectors<-subset.summary.dataframe[,grep("^q\\d",names(subset.summary.dataframe))]
		bestModelIndex=which(subset.summary.dataframe$deltaAIC==0)[1]

		meanBirthVector<-rep(NA,dim(allBirthVectors)[2])
		for (i in 1:length(meanBirthVector)) {
			meanBirthVector[i]=weighted.mean(x=allBirthVectors[,i],w=subset.summary.dataframe$AICweight,na.rm=TRUE)
		}
		names(meanBirthVector)<-paste("modelavg_",names(allBirthVectors),sep="")

		bestBirthVector<-as.numeric(allBirthVectors[bestModelIndex,])
		names(bestBirthVector)<-paste("best_",names(allBirthVectors),sep="")


		meanDeathVector<-rep(NA,dim(allDeathVectors)[2])
		for (i in 1:length(meanDeathVector)) {
			meanDeathVector[i]=weighted.mean(x=allDeathVectors[,i],w=subset.summary.dataframe$AICweight,na.rm=TRUE)
		}
		names(meanDeathVector)<-paste("modelavg_",names(allDeathVectors),sep="")
		
		bestDeathVector<-as.numeric(allDeathVectors[bestModelIndex,])
		names(bestDeathVector)<-paste("best_",names(allDeathVectors),sep="")


		bestDiversificationVector<-bestBirthVector-bestDeathVector
		meanDiversificationVector<-meanBirthVector-meanDeathVector
		for (charStateRow in 1:(2^nchar)) { 
			binaryStateRowVector<-digitsBase(charStateRow-1,ndigits=nchar)[,1]
			names(bestDiversificationVector)[charStateRow]<-paste("best_netdivers",paste(binaryStateRowVector,sep="",collapse=""),sep="")
			names(meanDiversificationVector)[charStateRow]<-paste("modelavg_netdivers",paste(binaryStateRowVector,sep="",collapse=""),sep="")
		}


		meanTransitionVector<-rep(NA,dim(allTransitionVectors)[2])
		for (i in 1:length(meanTransitionVector)) {
			meanTransitionVector[i]=weighted.mean(x=allTransitionVectors[,i],w=subset.summary.dataframe$AICweight,na.rm=TRUE)
		}
		names(meanTransitionVector)<-names(allTransitionVectors)
		
		bestTransitionVector<-as.numeric(allTransitionVectors[bestModelIndex,])
		names(bestTransitionVector)<-names(allTransitionVectors) #don't want "best" prefix here or matching later won't work (well, will need modification to work)


		meaningfulMeanTransitionVector<-rep(NA,length(qNamesVector))
		meaningfulBestTransitionVector<-rep(NA,length(qNamesVector))
		for (i in 1:length(qNamesVector)) {
			meaningfulMeanTransitionVector[i]<-meanTransitionVector[which(names(meanTransitionVector)==qNamesVector[i])]
			meaningfulBestTransitionVector[i]<-bestTransitionVector[which(names(bestTransitionVector)==qNamesVector[i])]
		}
		names(meaningfulMeanTransitionVector)<-paste("modelavg_",qNamesVector,sep="")
		names(meaningfulBestTransitionVector)<-paste("best_",qNamesVector,sep="")
		tmp.matrix=matrix(c(subset.summary.dataframe[bestModelIndex,1:7], subset.summary.dataframe[bestModelIndex,11:12], subset.summary.dataframe[bestModelIndex,"AICweight"],bestBirthVector,bestDeathVector,bestDiversificationVector,meaningfulBestTransitionVector,meanBirthVector,meanDeathVector,meanDiversificationVector,meaningfulMeanTransitionVector),ncol=1)
		row.names(tmp.matrix)=c(names(subset.summary.dataframe[bestModelIndex,1:7]),names(subset.summary.dataframe[bestModelIndex,11:12]),"AICweight",names(bestBirthVector),names(bestDeathVector),names(bestDiversificationVector),names(meaningfulBestTransitionVector),names(meanBirthVector),names(meanDeathVector),names(meanDiversificationVector),names(meaningfulMeanTransitionVector))
		if (dim(summary2.dataframe)[1]==0) {
			summary2.dataframe=data.frame(unlist(tmp.matrix),row.names=row.names(tmp.matrix))
		}
		else {
			summary2.dataframe<-cbind(summary2.dataframe,data.frame(unlist(tmp.matrix),row.names=row.names(tmp.matrix)))
		}
		names(summary2.dataframe)[dim(summary2.dataframe)[2] ]=partitionSchemeText
		write.table(summary2.dataframe,file="../Summaries/OverallSummary.txt",sep="\t")
		save(summary2.dataframe,file="../Summaries/OverallSummary.Rsave",ascii=TRUE)

	}
}
write.table(summary2.dataframe,file="../Summaries/OverallSummary.txt",sep="\t")
save(summary2.dataframe,file="../Summaries/OverallSummary.Rsave",ascii=TRUE)

