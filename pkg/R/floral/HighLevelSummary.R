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
		save(highlevel.dataframe,file="/Users/bomeara/Sites/Floral/RunsApril2011/Summaries/Highlevel.dataframe.Rsave",compress=TRUE)
		print(paste("just did  T",T,"D",D,"with length =",dim(highlevel.dataframe)[1]))
	}
}

ls()
deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
relativeLikelihood<-exp(-0.5 * deltaAIC)
AICweight<-relativeLikelihood/sum(relativeLikelihood)
highlevel.dataframe<-cbind(deltaAIC,AICweight,highlevel.dataframe)
save(highlevel.dataframe,file="/Users/bomeara/Sites/Floral/RunsApril2011/Summaries/Highlevel.dataframe.Rsave",compress=TRUE)
summarizeModelWeights(summary.dataframe=highlevel.dataframe,S=S,transitionModels=transitionModels, diversificationModels=diversificationModels)
