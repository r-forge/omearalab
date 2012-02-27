source("V6_UtilityFns.R")
setwd("../Summaries")
ls()
highlevel.dataframe<-data.frame()
load(file="/Users/bomeara/Sites/RunsJan2012/Summaries/Highlevel.dataframe.withrates.Rsave")
ls()
toDelete<-c()
for (i in 1:dim(highlevel.dataframe)[1]) {
	if(length(grep(".x..xx",highlevel.dataframe$focal[i],perl=TRUE))!=1) {
		toDelete<-c(toDelete,-1*i)
	}
}
highlevel.dataframe<-highlevel.dataframe[toDelete,]
deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
relativeLikelihood<-exp(-0.5 * deltaAIC)
AICweight<-relativeLikelihood/sum(relativeLikelihood)
highlevel.dataframe<-cbind(deltaAIC,AICweight,highlevel.dataframe)
save(highlevel.dataframe,file="/Users/bomeara/Sites/RunsJan2012/Summaries/Highlevel.dataframe.withrates.noBoring.Rsave",compress=TRUE)
summarizeModelWeights(summary.dataframe=highlevel.dataframe,S=S,transitionModels=transitionModels, diversificationModels=diversificationModels)
