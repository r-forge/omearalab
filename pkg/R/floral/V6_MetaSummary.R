setwd("/Users/bomeara/Sites/RunsJan2012/Summaries")
load("Highlevel.dataframe.withrates.Rsave")
highlevel.dataframe<-highlevel.dataframe[,c(-3,-4)] #dropping misleading columns
best<-which.min(highlevel.dataframe$AIC)
print(highlevel.dataframe[best,1:13])
plot(highlevel.dataframe$k_all,highlevel.dataframe$deltaAIC+1,log="y")
print(dim(highlevel.dataframe))
print("final dim should be 21126  4173")