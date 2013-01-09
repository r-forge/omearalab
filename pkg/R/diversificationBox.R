library(Rmpfr)
prec<-400
t<-160 #units of millions of years
min.N<-250000
max.N<-500000

r.conversion<-function(b, eps) {
  return(b*(1-eps)) 
}

beta<-function(r, eps, t) {
  #return(exp(r*t + log(1 - exp(-r*t)) - r*t - log(1 - eps * exp(-r*t))))
  r<-mpfr(r, prec)
  t<-mpfr(t, prec)
  #eps<-mpfr(eps, prec)
  return((exp(r*t)-1)/(exp(r*t) - eps)) 
}

alpha<-function(r, eps, t) {
  return(eps * beta(r, eps, t)) 
}



pN<-function(N, eps, r, t) {
  return(beta(r, eps, t)^(N-1) - beta(r, eps, t)^(N))
}

negloglikelihood<-function(x, b, N, t, bad.value=100000000) {
  r=r.conversion(b, as.numeric(x))
  prob<-pN(N, eps=x, r, t)
  neglogl<- -log(prob)
  #print(as.numeric(c(x, prob, neglogl)))
  if(!is.finite(neglogl)) {
    neglogl<-bad.value
  }
  return(neglogl)
}


n.species.vector<-seq(from=min.N, to = max.N, length.out=10)
b.vector<-seq(from=log(min(n.species.vector))/t, to=1, length.out=100)
results<-matrix(nrow=length(n.species.vector), ncol=length(b.vector))
for (n.species.index in sequence(length(n.species.vector))) {
  N<-n.species.vector[n.species.index]
  for (b.index in sequence(length(b.vector))) {
    b<-b.vector[b.index]
    result<-optimizeR(f=negloglikelihood,lower=0, upper=1.1, b=b, N=N, t=t)
    print(c(N, b, as.numeric(result$minimum)))
    results[n.species.index, b.index]<-as.numeric(result$minimum)
  }
}
rownames(results)<-paste(round(n.species.vector)," species")
colnames(results)<-paste(round(b.vector, 4)," birth")
print(results)
pdf(file="/Users/bomeara/Dropbox/CollabBeaulieu/angiosperms.div.pdf")
contour(x=n.species.vector, y=b.vector, z=results,xlab="n.species",ylab="speciation rate", bty="n")
dev.off()