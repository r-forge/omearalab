library(Rmpfr)
prec<-400
t<-50 #units of millions of years
min.N<-50
max.N<-50

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
b.vector<-seq(from=log(min(n.species.vector))/t, to=1, length.out=50)
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
pdf(file="~/Dropbox/CollabBeaulieu/metazoa.pdf")
contour(x=n.species.vector, y=b.vector, z=results,xlab="n.species",ylab="speciation rate absolute", bty="n", nlevels=25)
axis(side=4, at=seq(from=min(b.vector), to=max(b.vector), length.out=5), labels=round(seq(from=min(b.vector), to=max(b.vector), length.out=5)/(log(max(n.species.vector))/t),2))

dev.off()



###Eq.6 Magallon and Sanderson 2001:
n=50
e=.75
t=50
rc<-(1/t)*(log(.5*n*(1-e^2)+2*e+.5*(1-e)*sqrt(n*(n*(e^2)-8*e+2*n*e+n)))-log(2))
