library(ape)
library(phytools)
library(geiger)
library(picante)
library(parallel)
source("~/Documents/MyDocuments/Active/OMearaLabR/pkg/R/PhyloWizard.R")

true.phy <- read.tree("~/Dropbox/CollabBeaulieu/MakingUpData/Apiales_true.tre") #has about 1500 tips. 
true.phy$edge.length <- true.phy$edge.length / max(branching.times(true.phy))
taxonomy.phy <- read.tree("~/Dropbox/CollabBeaulieu/MakingUpData/Apiales_taxonomyTOTAL.tre")
nrep<-5
true.cloud<-list(true.phy)
taxonomy.cloud <- PhyloWizard(tips=taxonomy.phy$tips, constraint=taxonomy.phy, nrep=nrep)
model.names<-c("True","PhyloWizardTaxonomy")

all.clouds<-list(true.cloud, taxonomy.cloud)
print("Heights")
print(RunAnalysis(all.clouds, GetHeight))

sampling.vector <- c(.05, .25, .50, .75, .95)

#subsampling only, using pd weight to bias
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, SubsampleTaxa(true.phy, pd.weight=1, f=sampling.vector[i]), simplify=FALSE)))
  model.names<-append(model.names, paste("SubsampleOnly_PDWeight_", sampling.vector[i], sep=""))
}
print("Heights")
print(RunAnalysis(all.clouds, GetHeight))

#subsampling, uniform probability
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, SubsampleTaxa(true.phy, pd.weight=0, f=sampling.vector[i]), simplify=FALSE)))
  model.names<-append(model.names, paste("SubsampleOnly_UnifWeight_", sampling.vector[i], sep=""))
}
print("Heights")
print(RunAnalysis(all.clouds, GetHeight))

#subsampling, then using taxonomy to fill in missing, using pd weight to bias subsampling
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, DoSingleResolve(CombineTaxonomyAndSubsampledTrees(taxonomy.phy, SubsampleTaxa(true.phy, pd.weight=1, f=sampling.vector[i]))), simplify=FALSE)))
  model.names<-append(model.names, paste("SubsamplePlusTaxonomy_PDWeight_", sampling.vector[i], sep=""))
  
}
print("Heights")
print(RunAnalysis(all.clouds, GetHeight))

#subsampling, then using taxonomy to fill in missing, using unbiased subsampling
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  list(replicate(nrep, DoSingleResolve(CombineTaxonomyAndSubsampledTrees(taxonomy.phy, SubsampleTaxa(true.phy, pd.weight=0, f=sampling.vector[i]))), simplify=FALSE)))
  model.names<-append(model.names, paste("SubsamplePlusTaxonomy_UnifWeight_", sampling.vector[i], sep=""))
  
}
print("Heights")
print(RunAnalysis(all.clouds, GetHeight))

save(all.clouds, true.phy, taxonomy.phy, nrep, sampling.vector, file=paste("~/Dropbox/CollabBeaulieu/MakingUpData/nrep", nrep, ".RData", sep=""))

Simulate.BM <- function(phy) {
  return(fastBM(phy, sig2=1, a=1)) 
}

Simulate.Binary <- function(phy) {
  return(sim.char(phy, model.matrix=list(rbind(c(-.6, .6), c(.1, -.1))), model="discrete", n=100))
}

Analyze.BM <- function(phy, traits) {
  cleaned<-treedata(phy, traits)
  return(fitContinuous(cleaned$tree, cleaned$data, model="BM")$Trait1$beta)
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

GetHeight <- function(phy) {
  return(max(branching.times(phy))) 
}

Height.results<-RunAnalysis(all.clouds, GetHeight)
save(list=ls(), file=paste("~/Dropbox/CollabBeaulieu/MakingUpData/Height.nrep", nrep, ".RData", sep=""))
print(cbind(model.names, sapply(Height.results, mean, na.rm=TRUE)))


PD.results<-RunAnalysis(all.clouds, Analyze.PD.Sample)
save(list=ls(), file=paste("~/Dropbox/CollabBeaulieu/MakingUpData/PD.nrep", nrep, ".RData", sep=""))
print(cbind(model.names, sapply(Height.results, mean, na.rm=TRUE), sapply(PD.results, mean, na.rm=TRUE)))
            

traits<-Simulate.BM(true.phy)
BM.results<-RunAnalysis(all.clouds, Analyze.BM, traits=traits)
save(list=ls(), file=paste("~/Dropbox/CollabBeaulieu/MakingUpData/PDandBM.nrep", nrep, ".RData", sep=""))
print(cbind(model.names, sapply(Height.results, mean, na.rm=TRUE), sapply(PD.results, mean, na.rm=TRUE), sapply(BM.results, mean, na.rm=TRUE)))
