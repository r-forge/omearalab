
maxsteps<-100000
hittingTimes<-c()
nreps<-2000
for (i in sequence(nreps)) {
  state=rnorm(1,0,1)
  cross=FALSE
  step<-0
  while (!cross && step<maxsteps) {
    step=step+1
    newstate<-rnorm(1,state,1)
    if (sum(abs(state),abs(newstate))!=abs(state+newstate)) {
      cross=TRUE
      hittingTimes<-append(hittingTimes,step)
    }
    else {
      state<-newstate 
      if (step==maxsteps) {
        hittingTimes<-append(hittingTimes,NA)
      }
    }
  }
}
plot(density(hittingTimes,na.rm=TRUE))

nreps<-5000
timeToCross0<-function(maxsteps=100000,state=rnorm(1,0,1)) {
  cross=FALSE
  step<-0
  while (!cross && step<maxsteps) {
    step=step+1
    newstate<-rnorm(1,state,1)
    if (sum(abs(state),abs(newstate))!=abs(state+newstate)) {
      cross=TRUE
      return(step)
    }
    else {
      state<-newstate 
      if (step==maxsteps) {
        return(NA)
      }
    }
  }
}

times<-simplify2array(mclapply(rep(maxsteps,nreps),timeToCross0))
plot(density(times,na.rm=TRUE))
mean(times,na.rm=TRUE)

nreps<-5000
maxstepsBase<-10
numberDoublings<-17
meanTimes<-c()
medianTimes<-c()
sdTimes<-c()
maxstepVect<-c()
for (i in sequence(numberDoublings)) {
   times<-simplify2array(mclapply(rep(maxstepsBase*(2^i),nreps),timeToCross0))
   meanTimes<-append(meanTimes,mean(times,na.rm=TRUE))
   medianTimes<-append(medianTimes,median(times,na.rm=TRUE))
   sdTimes<-append(medianTimes,sd(times,na.rm=TRUE))
   maxstepVect<-append(maxstepVect, maxstepsBase*(2^i))
   print(c(i,meanTimes[i],medianTimes[i],sdTimes[i]))
}

plot(maxstepVect,meanTimes)


multipleHitTimes<-function(maxsteps=10000,state=rnorm(1,0,1)) {
  hitTimes<-c()
  totalSteps=0
  while(totalSteps<maxsteps) {
     hitTime<-timeToCross0(maxsteps-totalSteps)
     if (!is.na(hitTime)) {
       hitTimes<-append(hitTimes,hitTime+totalSteps)
      totalSteps<-totalSteps+hitTime

     }
     else {
        return(hitTimes) 
     }
  }
}

countMultipleHitTimes<-function(maxsteps=10000,state=rnorm(1,0,1)) {
  return(length(multipleHitTimes(maxsteps=maxsteps,state=state)))
}

maxsteps=100
times<-unlist(mclapply(rep(maxsteps, 100), multipleHitTimes,mc.cores=7))
plot(density(times,from=0,to=maxsteps))
plot(hist(times,breaks=c(-1:maxsteps+1)))
changes<-unlist(mclapply(rep(maxsteps, 100), countMultipleHitTimes,mc.cores=7))
mean(changes)
plot(hist(changes,breaks=c(-1:max(changes)+1)))