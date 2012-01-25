setwd("/Users/bomeara/Sites/RunsNov2011/Summaries")
load("Highlevel.dataframe.withrates.Rsave")
highlevel.dataframe<-highlevel.dataframe[,c(-3,-4)] #dropping misleading columns
best<-which.min(highlevel.dataframe$AIC)
print(highlevel.dataframe[best,1:13])
plot(highlevel.dataframe$k_all,highlevel.dataframe$AIC+1,log="y")
print(highlevel.dataframe[best,which(as.numeric(highlevel.dataframe[best,])>0.00000000000000001)])