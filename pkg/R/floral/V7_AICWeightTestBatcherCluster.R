setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
library(diversitree)
source("V6_UtilityFns.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_CommandsLocal.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/AICweightJan2013")
load(file="testing.model.subset.Rsave")
library(foreach)
library(doMC)
registerDoMC(13)



generating.models<-testing.model.subset[c(1, 2, 51, 51),]
generating.models[dim(generating.models)[1], sapply(generating.models[1,], is.numeric)]<-colMeans(testing.model.subset[, sapply(generating.models[1,], is.numeric)])
generating.models[dim(generating.models)[1], c(1, 2, 5, 7, 8, 9)]<-rep(NA, 6)
generating.models[dim(generating.models)[1], 3]<-"average"
nreps<-100

save(list=ls(), "AICWorkspace.RSave", compress=TRUE) #has all the local objects

foreach  (rep=sequence(nreps)) %dopar% try(doRunWithinForeach(rep, generating.models, testing.model.subset))

