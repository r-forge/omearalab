library(parallel)

source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V6_StochasticSSASims.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V6_NewSimulator.R")

files<-paste("simulation",sequence(500),".RSave",sep="")
mclapply(files, doParallelSSA, x0=x0, a=a, nu=nu, parms=parms, mc.cores=detectCores())
