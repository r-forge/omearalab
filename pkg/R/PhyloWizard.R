require(ape)
require(phytools)


#THIS IS INTENDED AS A SNARKY RESPONSE TO APPROACHES THAT ALLOW PEOPLE TO
#  RESPOND TO LACK OF INFORMATION ABOUT A TREE BY MAKING IT UP
#  IT ALLOWS YOU TO DO THIS, PERHAPS WITH A CONSTRAINT TREE
#  IN THIS CASE, THE CONSTRAINT TREE MUST HAVE ALL THE TAXA, BUT
#  BE INCOMPLETELY RESOLVED (THINK OF THE NCBI TAXONOMY TREE)
#  THE RESULTING TREE WILL HAVE REASONABLE-LOOKING BRANCH LENGTHS
#  IF THE CONSTRAINT TREE HAS A MEANINGFUL DEPTH; OTHERWISE, THE DEPTH
#  WILL SIMPLY BE ONE.

#  Brian O'Meara, Dec. 4, 2012

#example: 
# phy<-di2multi(rcoal(10), 0.5) # check to make sure this is incompletely resolved

# trees<-PhyloWizard(constraint=phy, nrep=10) 
# plot(trees[[1]]) #note that the resulting object is a list of phylo objects, not multiPhylo
# plot(trees[[2]])

#this comes from phytools.
#modified to use a taxonomy tree with all taxa as a starting tree
mrp.supertree.modified<-function (phy.full, phy.sub, n.best.observed=10) 
{
  phy<-c(phy.full, phy.sub)
  if (!class(phy) == "multiPhylo") 
    stop("phy must be object of class 'multiPhylo.'")
  X <- list()
  characters <- 0
  for (i in 1:length(phy)) {
    temp <- prop.part(phy[[i]])
    X[[i]] <- matrix(0, nrow = length(phy[[i]]$tip), ncol = length(temp) - 
                       1)
    for (j in 1:ncol(X[[i]])) X[[i]][c(temp[[j + 1]]), j] <- 1
    rownames(X[[i]]) <- attr(temp, "labels")
    if (i == 1) 
      species <- phy[[i]]$tip.label
    else species <- union(species, phy[[i]]$tip.label)
    characters <- characters + ncol(X[[i]])
  }
  XX <- matrix(data = "?", nrow = length(species), ncol = characters, 
               dimnames = list(species))
  j <- 1
  for (i in 1:length(X)) {
    XX[rownames(X[[i]]), c(j:((j - 1) + ncol(X[[i]])))] <- X[[i]][1:nrow(X[[i]]), 
                                                                  1:ncol(X[[i]])]
    j <- j + ncol(X[[i]])
  }
  contrast <- matrix(data = c(1, 0, 0, 1, 1, 1), 3, 2, dimnames = list(c("0", 
                                                                         "1", "?"), c("0", "1")), byrow = TRUE)
  XX <- phyDat(XX, type = "USER", contrast = contrast)
  best.score <- parsimony(phy.full, data=XX)
  best.tree <- multi2di(phy.full)
  n.best.seen <- 0
  while(n.best.seen < n.best.observed) { #keep searching until you get trees of the best observed score n.best.observed times
    print("doing pratchet")
    start.phy <- multi2di(phy.full)
    start.index <- runif(1,0,1)
    if (start.index<0.4) {
      start.phy <- best.tree 
      print("starting from best tree")
    } 
    if (start.index>0.4 & start.index < 0.8) {
      start.phy <- random.addition(XX)
      print("starting from random.addition tree")
    }
    
    supertree <- pratchet(XX, start=start.phy, trace = 0, all = FALSE)
    if (best.score == attr(supertree, "pscore")) {
      n.best.seen <- n.best.seen + 1
    }
    if (best.score > attr(supertree, "pscore")) {
      n.best.seen <- 0
      print(paste("replacing tree of score ", best.score, "with tree of score ", attr(supertree, "pscore")))
      best.score <- attr(supertree, "pscore")
      best.tree <- multi2di(di2multi(acctran(supertree, XX), 0.01)) 
    }
  }
  return(best.tree)
}

CombineTaxonomyAndSubsampledTrees <- function(taxonomy.phy, subtree.phy) {
  
  return(mrp.supertree.modified(taxonomy.phy, subtree.phy))
}

PhyloWizard <- function( tips=NA, constraint=NA, nrep=100, verbose=TRUE) {
  phy<-NULL
  if (verbose) {
    print("now making up branch lengths and trees") 
  }
  if(!is.na(constraint[1])) {
    
    return(replicate(nrep, DoSingleResolve(constraint, verbose), simplify=FALSE))
  }
  if (length(tips)>2) {
    return(replicate(nrep, DoSingleResolve(stree(n=length(tips), type="star", tip.label=tips), verbose), simplify=FALSE))
  }
  stop("Cannot make up trees until you specify either taxon names or constraint trees") 
}

DoSingleResolve <- function(constraint, verbose=TRUE) {
  if(verbose) {
    print("resolving tree")
    flush.console()
  }
  phy<-multi2di(constraint,  random=TRUE) #randomly resolves polytomies
  if(verbose) {
    print("making up branch lengths")
    flush.console()
  }
  return(MakeUpBranchlengths(phy, verbose))
}

MakeUpBranchlengths <- function(phy, verbose=TRUE) {
  if(verbose) {
    warning(paste("Responding to lack of information (in this case, lack of information about the phylogeny) by making it up is generally considered a bad thing"))
  }
  tree.depth<-1
  if (is.null(phy$edge.length)) {
    phy<-compute.brlen(phy, method=0)
  }
  else {
    tree.depth<-max(branching.times(phy)) 
  }
  expected.rate<-log(Ntip(phy)/2)/tree.depth #expected speciation rate under a Yule
  random.brlen<-rexp(length(which(phy$edge.length==0)), expected.rate) #draws brlen from exponential wait times
  phy$edge.length[which(phy$edge.length==0)]<-random.brlen 
  phy<-chronopl(phy, lambda=0, age.min=max(branching.times(phy)), age.max=max(branching.times(phy))) #makes the tree ultrametric, making as few changes as possible to the given brlen
  phy$edge.length <- phy$edge.length / max(branching.times(phy))
  return(phy)
}


#subsample f, the fraction of taxa to save, in a tree. 
#pd.weight=0 -> all taxa have equal chance of being sampled
#pd.weight=1 -> taxa sampled based on pd; those with more pd have higher chance of being sampled
#these weights can be varied: pd.weight=0.5, for example
SubsampleTaxa <- function(phy, pd.weight=0, f=0.5) {
  pd.vector <- phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  pd.vector <- pd.vector/sum(pd.vector)
  pd.vector <- 1-pd.vector
  tips.ordered <- phy$tip.label[phy$edge[which(phy$edge[,2]<=Ntip(phy)),2]]
  tips.doomed <- sample(tips.ordered, size=round((1-f)*Ntip(phy)), replace=FALSE, prob=pd.vector*pd.weight+(1-pd.weight)*rep(1/Ntip(phy), Ntip(phy)))
  return(drop.tip(phy, tips.doomed))
}

#p is the fraction of shortest edges you want to delete. A value of 0.4 will leave the tree with 60% of its original internal nodes
ShrinkByPercent<-function(phy, p) {
  internal.edge.length<-phy$edge.length[which(phy$edge[,2]>Ntip(phy))]
  cutoff<-quantile(internal.edge.length, probs=p)
  return(di2multi(phy, cutoff))
}