library(GillespieSSA)

load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/Summaries/Highlevel.dataframe.withrates.Rsave")
names(highlevel.dataframe)[1:10]
highlevel.dataframe$deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
highlevel.dataframe$AICweight<-exp( -0.5 * highlevel.dataframe$deltaAIC)
highlevel.dataframe$AICweight <- highlevel.dataframe$AICweight/sum(highlevel.dataframe$AICweight)
highlevel.dataframe<-highlevel.dataframe[order(highlevel.dataframe$AICweight, decreasing=TRUE),]

key.focal.vector<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")

getMatchingTransitions<-function(focalRow, key.focal.from, key.focal.to, change.positions=c(2, 5, 6), average=TRUE) {
  key.focal.from.split<-strsplit(key.focal.from, "")[[1]]
  key.focal.to.split<-strsplit(key.focal.to, "")[[1]]
  matching.values<-c()
  for (i in c(0:1)) {
    key.focal.from.split[change.positions[1]]<-i
    key.focal.to.split[change.positions[1]]<-i   
    for (j in c(0:1)) {
      key.focal.from.split[change.positions[2]]<-j
      key.focal.to.split[change.positions[2]]<-j
      for (k in c(0:1)) {
        key.focal.from.split[change.positions[3]]<-k
        key.focal.to.split[change.positions[3]]<-k
        #print(names(focalRow)[1:25])
        #print(paste("q", paste(key.focal.from.split, sep="", collapse=""), "_", paste(key.focal.to.split, sep="", collapse=""), sep=""))
        matching.values<-append(matching.values, 
              focalRow[ which(names(focalRow)==paste("q", paste(key.focal.from.split, sep="", collapse=""), "_", paste(key.focal.to.split, sep="", collapse=""), sep=""))])
      }
    }
  }
  if (!average) {
    return(matching.values)
  }
  else {
    return(mean(unlist(matching.values)))
  }
}

q.values<-matrix(nrow=dim(highlevel.dataframe)[1],ncol=24)
q.names<-c()
for (focal.index.1 in sequence(length(key.focal.vector))) {
  for (focal.index.2 in c(focal.index.1:length(key.focal.vector))) {
    if (focal.index.1 != focal.index.2) {
      if( vectorMismatch(unlist(strsplit(key.focal.vector[focal.index.1],split="")), unlist(strsplit(key.focal.vector[focal.index.2], split="")))==1) {
        q.names<-append(q.names,paste("q", key.focal.vector[focal.index.1], "_", key.focal.vector[focal.index.2], sep=""))
        #print(paste("q", key.focal.vector[focal.index.1], "_", key.focal.vector[focal.index.2], sep=""))
        for (element in sequence(dim(highlevel.dataframe)[1])) {
          q.values[element, length(q.names)]<-getMatchingTransitions(highlevel.dataframe[1,], key.focal.vector[focal.index.1], key.focal.vector[focal.index.2])
        }
       # print(getMatchingTransitions(highlevel.dataframe[1,], key.focal.vector[focal.index.1], key.focal.vector[focal.index.2]))
        #print(apply(highlevel.dataframe[1:10,], 1, getMatchingTransitions, key.focal.from = key.focal.vector[focal.index.1], key.focal.to = key.focal.vector[focal.index.2] ))
        #q.values<-cbind(q.values,data.frame(matrix(apply(highlevel.dataframe[1:10,], 1, getMatchingTransitions, key.focal.from = key.focal.vector[focal.index.1], key.focal.to = key.focal.vector[focal.index.2] ) , ncol=1)))
      }
    }
  }
}

#do matching like that below to get all the columns for q01001_01000 and so forth. Average the relevant columns for each row, 
#  then do weighted average of these across the rows based on AIC weights
#  for doing models with equal diversification, filter for those models that have the simplest diversification model and redo.

#grepl(key.focal.xtoperiod[2],names(highlevel.dataframe),perl=TRUE)

















printUniqueValues<-function(focalRow) {
  combo<-focalRow[which(names(focalRow)=="focal")]
  print(combo)
  focalRow.unnamed<-as.vector(focalRow)
  names(focalRow.unnamed)<-NULL
  combo.xto0<-gsub("x","0",combo)
  combo.mostxto0<-gsub("x","0",sub("x","1",combo)) #change first x to 1, others to 0
  combo.mismatch<-""
  if (length(grep("0",combo))>0) { #has a zero to move
    combo.mismatch<-sub("0","1",combo) #won't match combo
  }
  else { #has a 1 to move
    combo.mismatch<-sub("1","0",combo) #won't match combo
    
  }
  combo.mismatch.xto0<-gsub("x","0",combo.mismatch)
  combo.mismatch.mostxto0<-gsub("x","0",sub("x","1",combo.mismatch))
#  print(c(combo,combo.xto0, combo.mostxto0, combo.mismatch, combo.mismatch.xto0, combo.mismatch.mostxto0))
  #now get matches and mismatches: q_combo.xto0_combo.mismatch.xto0 is focal to nonfocal
  rates<-rep(NA,10)
  names(rates)<-c("lambdaF","muF","divF","lambdaN","muN","divN","qNF","qFN","qNN","qFF")
  try(rates[1]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("lambda",combo.xto0,sep=""))][[1]]))
  try(rates[2]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("mu",combo.xto0,sep=""))][[1]]))
  try(rates[3]<-rates[1]-rates[2])
  try(rates[4]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("lambda",combo.mismatch.xto0,sep=""))][[1]]))
  try(rates[5]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("mu",combo.mismatch.xto0,sep=""))][[1]]))
  try(rates[6]<-rates[4]-rates[5])
  try(rates[7]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("q",combo.mismatch.xto0,"_",combo.xto0,sep=""))][[1]]))
  try(rates[8]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("q",combo.xto0,"_",combo.mismatch.xto0,sep=""))][[1]]))
  try(rates[9]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("q",combo.mismatch.xto0,"_",combo.mismatch.mostxto0,sep=""))][[1]]))
  try(rates[10]<-as.numeric(focalRow.unnamed[which(names(focalRow)==paste("q",combo.xto0,"_",combo.mostxto0,sep=""))][[1]]))
  #print(rates)
  return(rates)
}


unique.rates<-t(apply(highlevel.dataframe,1,printUniqueValues))



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