#source("PhyloWizard.R")

library(geiger)



fractions<-seq(from=0, to=1, length.out=6)
ntax<-100
b<-1
d<-0.5
nrep<-10
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

par(mfcol=c(1,2))
plot(fractions, est.eps, type="l", bty="n", xlab="Fraction unknown", ylab="eps")
lines(x=range(fractions), y=rep(d/b, 2), lty="dotted")
plot(fractions, est.div, type="l", bty="n", xlab="Fraction unknown", ylab="div")
lines(x=range(fractions), y=rep(b-d, 2), lty="dotted")