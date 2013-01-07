source("PhyloWizard.R")

library(geiger)

#p is the fraction of shortest edges you want to delete. A value of 0.4 will leave the tree with 60% of its original internal nodes
ShrinkByPercent<-function(phy, p) {
   internal.edge.length<-phy$edge.length[which(phy$edge[,2]>Ntip(phy))]
   cutoff<-quantile(internal.edge.length, probs=p)
   return(di2multi(phy, cutoff))
}

fractions<-seq(from=0, to=1, length.out=5)
ntax<-100
b<-1
d<-0.5
nrep<-4
est.eps<-rep(NA,length(fractions))
est.div<-est.eps

for (fraction.index in sequence(length(fractions))) {
  eps.this.rep<-rep(NA, nrep)
  div.this.rep<-rep(NA, nrep)
  fraction<-fractions[fraction.index]
  for (rep in sequence(nrep)) {
    print(rep)
    print(fraction)
    result<-NA
    while(is.na(result)) {
      try(result<-birthdeath(PhyloWizard(constraint=ShrinkByPercent(prune.extinct.taxa(birthdeath.tree(b=b, d=d, taxa.stop=ntax, return.all.extinct=FALSE)), fraction), nrep=1)[[1]]))
    }
    eps.this.rep[rep]<-result$para[1]
    div.this.rep[rep]<-result$para[2]
  }
  est.eps[fraction.index]<-median(eps.this.rep)
  est.div[fraction.index]<-median(div.this.rep)
  
}

