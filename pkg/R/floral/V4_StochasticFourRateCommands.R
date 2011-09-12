library(GillespieSSA)



oneSim<-function(maxWallTime=Inf) {
 taskid<-Sys.getenv("SGE_TASK_ID")
  starttime<-Sys.time()
  x0<-c(state00=2,state01=0,state11=0,state10=0)

  qNN=0.00143651
  qNF=0.00600116
  qFN=0.00143651
  qFF=0.00143651
  birthN=5.893707
  deathN=5.856925
  birthF=2.03713
  deathF=2.0089
  
  #00=ancestor, 01=focal == xxx0xx1x
  parms<-c(
   q00_01=qNF,
   q00_10=qNN,
   q01_11=qFN,
   q01_00=qFN,
   q11_10=qNN,
   q11_01=qNF,
   q10_11=qNN,
   q10_00=qNN,
   birth00=birthN,
   birth01=birthF,
   birth11=birthN,
   birth10=birthN,
   death00=deathN,
   death01=deathF,
   death11=deathN,
   death10=deathN
  )
 


  aNames<-c(
    "00 to 01",
    "00 to 10",
    "01 to 00",
    "01 to 11",
    "11 to 01",
    "11 to 10",
    "10 to 11",
    "10 to 00",
    "birth00",
    "birth01",
    "birth11",
    "birth10",
    "death00",
    "death01",
    "death11",
    "death10"
  )
  
  a<-c(
      "state00*q00_01",
      "state00*q00_10",
      "state01*q01_00",
      "state01*q01_11",
      "state11*q11_01",
      "state11*q11_10",
      "state10*q10_11",
      "state10*q10_00",
      "state00*birth00",
      "state01*birth01",
      "state11*birth11",
      "state10*birth10",
      "state00*death00",
      "state01*death01",
      "state11*death11",
      "state10*death10"
  )
  
 #state order is 
#00
#01
#11
#10
  nu<-matrix(c(
    -1,1,0,0,
    -1,0,1,0,
    1,-1,0,0,
    0,-1,1,0,
    0,1,-1,0,
    0,0,-1,1,
    0,0,1,-1,
    1,0,0,-1,
    1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1,
    -1,0,0,0,
    0,-1,0,0,
    0,0,-1,0,
    0,0,0,-1
    ), byrow=FALSE,nrow=length(x0))
    
  #a has reactions: the possible step moves
 # aNames<-   c("F to N","N to F","F birth", "F death", "N birth", "N death")
#  a<-        c("F*qFN", "N*qNF", "F*birthF","F*deathF","N*birthN","N*deathN")
  #nu<-matrix(c(-1,1,     1,-1,    1,0,       -1,0,      0,1,       0,-1),byrow=FALSE,nrow=length(x0)) #rows are states (F,N) and cols are the processes

  out<-ssa(x0,a,nu,parms,tf=130,maxWallTime=maxWallTime,method="OTL")
  outvector<-c(as.numeric(Sys.time()-starttime,units="secs"),out$data[dim(out$data)[1],])
  names(outvector)[1]<-"CPUtime"
  names(outvector)[2]<-"Duration"
  cat(outvector,"\n",file=paste("StochasticFourRate_",taskid,".txt",sep=""),append=TRUE)
 # ssa.plot(out)
  return(outvector)
}

#IMPORTANT:
#The way the batching happens, each set of runs will end with a long run. So to estimate
#the proportion of sims that go extinct, you'll have to take this into account (something like
#using the wait time until a success to estimate the probability, rather than just taking the
#ratio of number of runs that ended up with species to the total number of runs). Another way 
#would be to calculate the raw ratio, but using just the first X simulations in each run, where
#X is the minimum number of simulations in any run


taskid<-Sys.getenv("SGE_TASK_ID")
overallstarttime<-Sys.time()
sims<-cbind(oneSim())
save(sims,file=paste("StochasticFourRate_",taskid,".RSave",sep=""),compress=TRUE)
minduration<-3600 #one hour
while(as.numeric(Sys.time()-overallstarttime,units="secs")<minduration) {
  sims2<- cbind(oneSim())
  if(!is.null(dim(sims2))) {
     sims<-cbind(sims,sims2)
     save(sims,file=paste("StochasticFourRate_",taskid,".RSave",sep=""),compress=TRUE)
  }
}