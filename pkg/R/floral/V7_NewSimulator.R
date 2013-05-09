library(RColorBrewer)

#tf = duration of the simulation. 
#full.history=TRUE stores the history of number of species in each category, but this can get HUGE with long sims. If turned off, just returns the final vector
#rescale.species, if not NULL, adjusts the birth and death rates to get, on average, the desired number of species
#    This is useful if you do a subsample of species in bisse to estimate rates but you want to simulate the full number of species
#yule.scale is a way of making it a pure birth model, a birth-death model, or something in between
# essentially, new lambda = lambda - yule.scale*mu and new mu = (1 - yule.scale) * mu
#    =0 means the input rates are used
#    =1 means that lambda is set to lambda - mu and mu is set to zero
#    other values in between 
#x0.rescale lets you decide what character frequency to use for scaling. If the root state has a low diversification rate, and you expect to
#    quickly move to a different set of states with higher div. rates, you may want to use these to rescale so you don't overshoot

#note history now has time elapsed
OMearaSSA<-function(x0, q.vector, lambda.vector, mu.vector, tf, maxWallTime, verbose=TRUE, print.freq=100, full.history=TRUE, rescale.species=NULL, yule.scale=0, history.steps.to.save=seq(from=1,to=floor(tf),length.out=floor(tf)), t.rescale=tf, x0.rescale=NULL) { #turn off full.history to save on memory
  #make sure order of states is the same in all input objects

  lambda.vector<-lambda.vector-(yule.scale*mu.vector)
  mu.vector<-mu.vector*(1-yule.scale)
  
  q.matrix<-matrix(0,nrow=length(x0),ncol=length(x0))
  rownames(q.matrix)<-names(x0)
  colnames(q.matrix)<-names(x0)
  time.elapsed<-0
  step.number<-0
  for (i in sequence(length(q.vector))) {
    q.matrix[ which(rownames(q.matrix)==extractFrom(names(q.vector)[i])), 
              which(colnames(q.matrix)==extractTo(names(q.vector)[i]))
              ]<-q.vector[i]
  }
  
  if(!is.null(rescale.species)) {
     #average.diversification<-weighted.mean(x = lambda.vector - mu.vector, w = x0) #average net diversification rate for starting taxa
     #ideal.rate<- (log( rescale.species / sum(x0) )) / t.rescale #this rate will give the desired number of species over the desired time
     #scale.factor<-ideal.rate / average.diversification
     if(is.null(x0.rescale)) {
     	x0.rescale<-x0
     }
    scale.factor<-getOptimalScaling(lambda=weightedHarmonicMeanZeroCorrection(x=lambda.vector, w=x0.rescale), mu=weightedHarmonicMeanZeroCorrection(x=mu.vector, w=x0.rescale), t=t.rescale, N0=sum(x0), N=rescale.species)
    print(paste("scale.factor is ",scale.factor))
    print("x0.rescale is ")
    print(x0.rescale)
    print("lambda") 
    print(cbind(lambda.vector, lambda.vector*scale.factor))
    print("mu")
    print(cbind(mu.vector, mu.vector*scale.factor))
     lambda.vector<-scale.factor * lambda.vector
     mu.vector<-scale.factor * mu.vector
  }
  
  if(verbose) {
    print(paste("The probability of extinction is roughly",probExtinction(weightedHarmonicMeanZeroCorrection(lambda.vector, x0), weightedHarmonicMeanZeroCorrection(mu.vector, x0), tf, sum(x0)))) 
  }
  
  names(lambda.vector)<-gsub("lambda","",names(lambda.vector))
  names(mu.vector)<-gsub("mu","",names(mu.vector))
  #print(q.matrix)
  all.matrix<-cbind(q.matrix,matrix(lambda.vector,ncol=1),matrix(mu.vector,ncol=1))
  if(verbose) {
    print(all.matrix)
  }
  start.time<-Sys.time()
  x<-x0
  history<-matrix(c(0,x0),nrow=1)
  history.step<-Inf
  history.index<-1
  if(length(history.steps.to.save)>0) {
  	history.step<-history.steps.to.save[history.index]
  }
  while ((Sys.time()-start.time)<maxWallTime) {
    step.number<-step.number+1
    rate.matrix<-ScaleTransitionProbs(all.matrix,x)
    time.interval<-rexp(1, sum(rate.matrix))
    
    if ((time.elapsed+time.interval)>tf) {
      history<-rbind(history,c((tf),x)) #done running
      return(history)
    }
    
    if ((time.elapsed+time.interval)>history.step) {
      history<-rbind(history,c((history.step),x))
      history.index<-history.index+1
      if(history.index>length(history.steps.to.save)) {
      	history.step<-Inf
      }
      else {
      	history.step<-history.steps.to.save[history.index]
      }
    }
    
    
    change.row<-which(rmultinom(1,1,apply(rate.matrix, 1, sum))==1)
    change.col<-which(rmultinom(1,1,rate.matrix[change.row,])==1)
    move<-""
    if(change.col<=length(x)) { #transition
      x[change.row]<-x[change.row]-1
      x[change.col]<-x[change.col]+1
      if (verbose) {
         move<-paste("t_",names(x0[change.row]),"->",names(x0[change.col]),sep="")
      }
    }
    if(change.col==(length(x)+1)) { #birth
      x[change.row]<-x[change.row]+1
      if (verbose) {
        move<-paste("b_",names(x0[change.row]),sep="")
      }
    }
    if(change.col==(length(x)+2)) { #death
      x[change.row]<-x[change.row]-1
      if (verbose) {
        move<-paste("d_",names(x0[change.row]),sep="")
      }
    }
   time.elapsed<-time.interval+time.elapsed
   if(full.history) {
      history<-rbind(history,c(time.elapsed,x))
   }
    if(verbose || ((step.number %% print.freq) ==0)) {
     print(paste(sum(x),paste(x,sep=" ",collapse=" "),time.interval, time.elapsed, move, sep=" "))
    }
    if(sum(x)==0) {
      if(verbose) {
        print("all extinct")
      }
     history<-rbind(history,c(time.elapsed,x))
     return(history) 
    }
  }
  return(history)
}

ScaleTransitionProbs<-function(all.matrix, x) {
  return(all.matrix*x) #this multiplies each element in the i-th row of the matrix by the i-th element of x.  
}

plot.OMearaSSA<-function(history) {
  ncombos<-dim(history)[2]-1
  mypalette<-brewer.pal(ncombos,"Dark2")
# plot(range(cumsum(history[,1])),range(history[,2:dim(history)[2] ]),xlab="time",ylab="ntaxa",type="n",bty="n") 
 plot(range((history[,1])),range(history[,2:dim(history)[2] ]),xlab="time",ylab="ntaxa",type="n",bty="n") #now has elapsed time 
 for (i in 2:dim(history)[2]) {
  #lines( cumsum(history[,1]), history[, i],col=mypalette[i-1])
  lines( (history[,1]), history[, i],col=mypalette[i-1]) #now has elapsed time
 }
  legend(x="topleft",legend=colnames(history)[2:(ncombos+1)], text.col=mypalette)
}

plot.histories<-function(histories) { #uses the final outcome, not the path along the way
  transformed.history<-t(histories[,2:9])
  ncombos<-dim(histories)[2]-1
  mypalette<-brewer.pal(ncombos,"Dark2")
  par(mfcol=c(2,1))
  barplot(transformed.history[, order(apply(transformed.history, 2, sum))],col=mypalette,border=NA, main="total species",space=0)
  legend(x="topleft",legend=colnames(histories)[2:(ncombos+1)], text.col=mypalette,cex=0.5,ncol=2)
  barplot(t(t(transformed.history)/(apply(transformed.history, 2, sum)))[, order(apply(transformed.history, 2, sum))],col=mypalette,border=NA,main="proportion",space=0) 
}

doParallelSSA<-function(x0, q.means, lambda.means, mu.means, tf=136, maxWallTime=Inf, verbose=F, file.string="", full.history=FALSE, print.freq=100, rescale.species=NULL, yule.scale=0, t.rescale=136, x0.rescale=NULL) {
  file.name<-paste("SSA_",file.string,"_",format(Sys.time(), "%b%d_%H_%M_%S"),"_",round(runif(1,1,1000000)),".RSave",sep="")
  survivors<-0
  history<-0
  attempts<-0
  while(survivors<=0) {
    attempts<-attempts+1
    history<-OMearaSSA(x0, q.means, lambda.means, mu.means, tf=tf, maxWallTime=maxWallTime, verbose=verbose,full.history=full.history, print.freq=print.freq, rescale.species=rescale.species, yule.scale=yule.scale, t.rescale=t.rescale, x0.rescale=x0.rescale)
    survivors<-history[dim(history)[1],2]
    print(c(attempts,max(apply(history[,2:9],1,sum))))
  }
  save(list=ls(), file=file.name)
}

appendParallelSSA<-function(x0, q.means, lambda.means, mu.means, prev.history, t.additional=1, maxWallTime=Inf, verbose=F, file.string="", full.history=FALSE, print.freq=100, rescale.species=NULL, yule.scale=0, t.rescale=136, x0.rescale=NULL) {
  file.name<-paste(file.string,"_","app",t.additional,"total", t.additional+max(prev.history[,1]), ".RSave",sep="")
  survivors<-0
  history<-0
  attempts<-0
  while(survivors<=0) {
    attempts<-attempts+1
    history<-OMearaSSA(x0=prev.history[dim(prev.history)[1],2:9], q.means, lambda.means, mu.means, tf=t.additional, maxWallTime=maxWallTime, verbose=verbose,full.history=full.history, print.freq=print.freq, rescale.species=rescale.species, yule.scale=yule.scale, t.rescale=t.rescale, x0.rescale=x0.rescale)
    survivors<-history[dim(history)[1],2]
    history<-rbind(prev.history,history[2:dim(history)[1],]) #the first row of the new history will be the same as the last of the old history
    history[,1]<-c(0:(dim(history)[1]-1)) #now get the right dates. Note this assumes 1 my spacing
    print(c(attempts,max(apply(history[,2:9],1,sum))))
  }
  save(list=ls(), file=file.name)
}


probExtinction<-function(lambda, mu, t, N0=2) {
  #after magallon and sanderson
  epsilon<-mu/lambda
  r<-lambda-mu
  beta<-(exp(r*t) - 1) / (exp(r*t) - epsilon)
  alpha<-epsilon * beta
  return(alpha^N0)
}


probGreaterThanEqualToNTaxaGivenObserved<-function(lambda, mu, t, N0=2, N=250000) {
  #after magallon and sanderson 2001, CORRECTED equation 11a
  epsilon<-mu/lambda
  r<-lambda-mu
  beta<-(exp(r*t) - 1) / (exp(r*t) - epsilon)
  alpha<-epsilon * beta
  k<-N
  result <- ((beta^(k-2)) / (1+alpha)) * ( k * (1 - alpha - beta + (alpha * beta) ) + alpha + 2*beta -1 )
  return(result)
}

probExactlyNTaxaNotConditioningOnObserved<-function(lambda, mu, t, N0=2, N=250000) {
  #after Raup 1985, equation A18, which is Bailey 1964, equation 8.47
  epsilon<-mu/lambda
  r<-lambda-mu
  beta<-(exp(r*t) - 1) / (exp(r*t) - epsilon)
  alpha<-epsilon * beta
  k<-N
  a<-N0
  individual.term<-function(x, a, alpha, beta, N) {
    return( choose(a, x) * choose(a+N-x-1, a-1) * (alpha^(a-x)) * (beta^(N-x)) * ((1 - alpha - beta)^x))
  }
  return(sum(sapply(c(0:min(a,N)), individual.term, a=a, alpha=alpha, beta=beta, N=N)))
}

probExactlyNTaxaGivenObserved<-function(lambda, mu, t, N0=2, N=250000) {
  epsilon<-mu/lambda
  r<-lambda-mu
  beta<-(exp(r*t) - 1) / (exp(r*t) - epsilon)
  alpha<-epsilon * beta
  a<-N0
  result<-probExactlyNTaxaNotConditioningOnObserved(lambda, mu, t, N0, N)/(1 - alpha^a)
  if(is.na(result)) {
    result<-0 
  }
  return(result)
}

likelihoodGivenScaling<-function(x, lambda, mu, t, N0=2, N=250000, multiplier=1) {
	return(multiplier*probExactlyNTaxaGivenObserved(x*lambda, x*mu, t, N0, N))
}

getOptimalScaling<-function(lambda, mu, t, N0=2, N=250000) {
  r<-lambda - mu
  ideal.rate<-log(N)/t
  scale.factor<-ideal.rate/r
  result<-optim(par=scale.factor, fn=likelihoodGivenScaling, method="L-BFGS-B", lambda=lambda, mu=mu, t=t, N0=N0, N=N, multiplier=-1)
  print(paste("Yule scale factor would be = ", scale.factor))
	return(result$par)
}