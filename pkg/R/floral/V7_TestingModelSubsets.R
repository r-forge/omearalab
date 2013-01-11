source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V6_UtilityFns.R")
load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/all.results.cleaned.Rsave")
testing.model.subset<-highlevel.dataframe[c(1:50, round(seq(from=100, to=19600, length.out=50))),] #100 models used for evaluating reproducibility, AIC weights
testing.model.subset$AICweight<-testing.model.subset$AICweight/sum(testing.model.subset$AICweight)
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/AICweightJan2013")
save(testing.model.subset, file="testing.model.subset.Rsave", compress=TRUE)
