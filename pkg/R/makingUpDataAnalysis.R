library(ape)
library(phytools)
library(geiger)

true.phy <- read.tree("~/Dropbox/CollabBeaulieu/MakingUpData/Apiales_true.tre") #has about 1500 tips. 
#taxonomy.phy <- read.tree("~/Dropbox/CollabBeaulieu/MakingUpData/Apiales_taxonomyTOTAL.tre")
#genera.phy <- Jeremy(taxonomy.phy)
#family.phy <- Jeremy(taxonomy.phy)
nreps<-5
true.cloud<-list(true.phy)
#taxonomy.cloud <- PhyloWizard(tips=taxonomy.phy$tips, constraints=taxonomy.phy, nreps=nreps)
#genera.cloud <- PhyloWizard(tips=taxonomy.phy$tips, constraints=genera.phy, nreps=nreps)
#family.cloud <- PhyloWizard(tips=taxonomy.phy$tips, constraints=family.phy, nreps=nreps)

#all.clouds<-list(true.cloud, taxonomy.cloud, genera.cloud, family.cloud)
all.clouds<-list(true.cloud)

sampling.vector <- c(.05, .25, .50, .75, .95)

#subsampling only, using pd weight to bias
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  replicate(nreps, SubsampleTaxa(true.phy, pd.weight=1, f=sampling.vector[i]), simplify=FALSE))
}

#subsampling, uniform probability
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  replicate(nreps, SubsampleTaxa(true.phy, pd.weight=0, f=sampling.vector[i]), simplify=FALSE))
}

#subsampling, then using taxonomy to fill in missing, using pd weight to bias subsampling
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  replicate(nreps, DoSingleResolve(CombineTaxonomyAndSubsampledTrees(c(taxonomy.phy, SubsampleTaxa(true.phy, pd.weight=1, f=sampling.vector[i])))), simplify=FALSE))
}

#subsampling, then using taxonomy to fill in missing, using unbiased subsampling
for (i in sequence(length(sampling.vector))) {
  print(sampling.vector[i])
  all.clouds<-append(all.clouds,  replicate(nreps, DoSingleResolve(CombineTaxonomyAndSubsampledTrees(c(taxonomy.phy, SubsampleTaxa(true.phy, pd.weight=0, f=sampling.vector[i])))), simplify=FALSE))
}

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

