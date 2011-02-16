library(diversitree) #obvious
library(sfsmisc) #for counting in binary
nchar=7
plot(x=c(0,1000000),y=c(0,8),type="n")

numberOfTransitionModelsForPartitionSize<-c(4,17,5+3*(2^3),5+3*(2^4),5+3*(2^5),5+3*(2^6),5+3*(2^7)) #make sure to include models that reduce to independent models: though they are present in smaller model sets, we might want independent transitions but dependent diversification
numberOfDiversificationModelsForPartitionSize<-c(19,12,6+3*(2^3),6+3*(2^4),6+3*(2^5),6+3*(2^6),6+3*(2^7))


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
for (partitionIterator in 1:((2^nchar)-1)) { #do not want a model with no chars, so start with 1, not zero
	#print(digitsBase(partitionIterator,ndigits=nchar))
	partitionScheme<-digitsBase(partitionIterator,ndigits=nchar)
	partitionScheme<-partitionScheme[,1] #1 means include that char in partition, NA means do not
	partitionScheme<-partitionScheme*c(1:nchar) #so this will result in a vector with digits showing if the character is included
	partitionSize<-length(which(partitionScheme>0))
	partitionSchemeText=paste(partitionScheme,sep="",collapse="_")
	for (transitionIndex in 1:numberOfTransitionModelsForPartitionSize[partitionSize]) { 
		for (diversificationIndex in 1:numberOfDiversificationModelsForPartitionSize[partitionSize]) {
			if (min(numberOfTransitionModelsForPartitionSize[partitionSize],numberOfDiversificationModelsForPartitionSize[partitionSize])>0) { #just to make sure we have models: for (1:1) and for (1:0) still run
				numberOfModels=numberOfModels+1
				#print(paste(partitionSize,numberOfModels))
				points(x=numberOfModels,y=partitionSize,pch=".",col="red")
			}
		}
	}
}
print(numberOfModels)

