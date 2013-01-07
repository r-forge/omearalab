require(ape)


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

PhyloWizard <- function( tips=NA, constraint=NA, nrep=100) {
  phy<-NULL
  if(!is.na(constraint[1])) {
    return(replicate(nrep, DoSingleResolve(constraint), simplify=FALSE))
  }
  if (length(tips)>2) {
    return(replicate(nrep, DoSingleResolve(stree(n=length(tips), type="star", tip.label=tips)), simplify=FALSE))
  }
  stop("Cannot make up trees until you specify either taxon names or constraint trees") 
}

DoSingleResolve <- function(constraint) {
  phy<-multi2di(constraint,  random=TRUE) #randomly resolves polytomies
  return(MakeUpBranchlengths(phy))
}

MakeUpBranchlengths <- function(phy) {
  warning(paste("Responding to lack of information (in this case, lack of information about the phylogeny) by making it up is generally considered a bad thing"))
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
   return(phy)
}