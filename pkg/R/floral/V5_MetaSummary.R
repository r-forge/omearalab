setwd("/Users/bomeara/Sites/RunsNov2011/Summaries")
load("Highlevel.dataframe.withrates.Rsave")
best<-which.min(highlevel.dataframe$AIC)
print(highlevel.dataframe[best,1:15])