library(ape)
library(phytools)
library(geiger)
library(picante)
library(parallel)
source("~/Documents/MyDocuments/Active/OMearaLabR/pkg/R/PhyloWizard.R")

true.phy <- read.tree("~/Dropbox/CollabBeaulieu/MakingUpData/Apiales_true.tre") #has about 1500 tips. 
taxonomy.phy <- read.tree("~/Dropbox/CollabBeaulieu/MakingUpData/Apiales_taxonomyTOTAL.tre")
nrep<-5
true.cloud<-list(true.phy)
taxonomy.cloud <- PhyloWizard(tips=taxonomy.phy$tips, constraint=taxonomy.phy, nrep=nrep)

all.clouds<-list(true.cloud, taxonomy.cloud)


sampling.vector <- c(.05, .25, .50, .75, .95)

#subsampling only, using pd weight to bias
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, SubsampleTaxa(true.phy, pd.weight=1, f=sampling.vector[i]), simplify=FALSE)))
}

#subsampling, uniform probability
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, SubsampleTaxa(true.phy, pd.weight=0, f=sampling.vector[i]), simplify=FALSE)))
}

#subsampling, then using taxonomy to fill in missing, using pd weight to bias subsampling
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, DoSingleResolve(CombineTaxonomyAndSubsampledTrees(taxonomy.phy, SubsampleTaxa(true.phy, pd.weight=1, f=sampling.vector[i]))), simplify=FALSE)))
}

#subsampling, then using taxonomy to fill in missing, using unbiased subsampling
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, DoSingleResolve(CombineTaxonomyAndSubsampledTrees(taxonomy.phy, SubsampleTaxa(true.phy, pd.weight=0, f=sampling.vector[i]))), simplify=FALSE)))
}

save(all.clouds, true.phy, taxonomy.phy, nrep, sampling.vector, file=paste("~/Dropbox/CollabBeaulieu/MakingUpData/nrep", nrep, ".RData", sep=""))

Simulate.BM <- function(phy) {
  return(sim.char(phy, model.matrix=matrix(1), nsims=1, model="brownian", root.state=1)) 
}

Simulate.Binary <- function(phy) {
  return(sim.char(phy, model.matrix=list(rbind(c(-.6, .6), c(.1, -.1))), model="discrete", n=100))
}

Analyze.BM <- function(phy, traits) {
  return(fitContinuous(phy, traits, model="BM")$Trait1$beta)
}

Analyze.Binary <- function(phy, traits) {
  return(fitDiscrete(phy, traits, model="ARD")$Trait1$q)
}

Analyze.BD <- function(phy) {
  return(birthdeath(phy)$para)
}

Analyze.PD.sample <- function(phy, nsamp=50) { #take a subsample of taxa on the tree, as you would to get a null PD. That way, trees with fewer taxa still have same expectation if they come from same overall tree.
  phy.sampled<-drop.tip(phy, sample(phy$tip.label, size = (Ntip(phy) - nsamp), replace=FALSE))
  return(sum(phy.sampled$edge.length))
}


RunAnalysis <- function(all.clouds, fn, ...) {
  all.results <- list()
  for (cloud in sequence(length(all.clouds))) {
    focal.cloud <- all.clouds[[cloud]]
    print(length(focal.cloud))
    ntrees <- length(focal.cloud)
    result <- rep(NA, ntrees)
    for (i in sequence(ntrees)) {
      if(class(all.clouds[[cloud]][[i]]) == "phylo") {
        result[i]<-fn(all.clouds[[cloud]][[i]], ...) 
      }
    }
    all.results<-append(all.results, list(result))
  }
  return(all.results)
}