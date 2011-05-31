source("V4_UtilityFns.R")
setwd("/Users/bomeara/Sites/Floral/RunsApril2011/Summaries")
ls()
highlevel.dataframe<-data.frame()
for (T in 1:5) {
	for (D in 1:6) {
		load(paste("RateSummaryT",T,"D",D,".Rsave"))
		if (T+D==2) { #first one
			highlevel.dataframe<-summary.dataframe[,3:13]
		}
		else {
			highlevel.dataframe<-rbind(highlevel.dataframe,summary.dataframe[,3:13])
		}
		save(highlevel.dataframe,file="Highlevel.dataframe.Rsave",compress=TRUE)
		print(paste("just did  T",T,"D",D,"with length =",dim(highlevel.dataframe)[1]))
	}
}

highlevel.dataframe->summary.dataframe
ls()
deltaAIC<-summary.dataframe$AIC-min(summary.dataframe$AIC)
relativeLikelihood<-exp(-0.5 * deltaAIC)
AICweight<-relativeLikelihood/sum(relativeLikelihood)
summary.dataframe<-cbind(deltaAIC,AICweight,summary.dataframe)
summarizeModelWeights(summary.dataframe=summary.dataframe,S=S,transitionModels=transitionModels, diversificationModels=diversificationModels)
