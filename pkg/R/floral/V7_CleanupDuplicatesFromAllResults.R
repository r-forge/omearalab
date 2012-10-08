setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/UnifiedApproachScripts")
load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/UnifiedApproachScripts/allresults.Rsave")

countX<-function(x) {
  return(sum(grepl("x",strsplit(x,"")[[1]])))
}

x.counts<-sapply(highlevel.dataframe$focal, countX)
highlevel.dataframe<-cbind(highlevel.dataframe, x.counts=x.counts)





CleanDuplicateModels<-function(highlevel.dataframe) {
  cleaned.highlevel.dataframe<-data.frame()
  possible.symmetry.T<-c(1,5)
  possible.symmetry.D<-c(1,3,6)
  impossible.symmetry.T<-c(2,3,4)
  impossible.symmetry.D<-c(2,4,5)
  for( T.index in sequence(length(impossible.symmetry.T))) {
    T<-impossible.symmetry.T[T.index]
    print(paste("T",T))
    cleaned.highlevel.dataframe<-rbind(cleaned.highlevel.dataframe, highlevel.dataframe[which(highlevel.dataframe$T==T),])
    highlevel.dataframe<-highlevel.dataframe[-which(highlevel.dataframe$T==T),]
  }
  for( D.index in sequence(length(impossible.symmetry.D))) {
    D<-impossible.symmetry.D[D.index]
    print(paste("D", D))
    
    cleaned.highlevel.dataframe<-rbind(cleaned.highlevel.dataframe, highlevel.dataframe[which(highlevel.dataframe$D==D),])
    highlevel.dataframe<-highlevel.dataframe[-which(highlevel.dataframe$D==D),]
  }
  #there is only one model with equal trans + yule, as focal combo doesn't matter
  cleaned.highlevel.dataframe<-rbind(cleaned.highlevel.dataframe, (highlevel.dataframe[which(highlevel.dataframe$D==1 & highlevel.dataframe$T==1),])[1,])
  highlevel.dataframe<-highlevel.dataframe[-which(highlevel.dataframe$D==1 & highlevel.dataframe$T==1),]
  
  #there is only one model with equal trans + birthdeath, as focal combo doesn't matter
  cleaned.highlevel.dataframe<-rbind(cleaned.highlevel.dataframe, (highlevel.dataframe[which(highlevel.dataframe$D==3 & highlevel.dataframe$T==1),])[1,])
  highlevel.dataframe<-highlevel.dataframe[-which(highlevel.dataframe$D==3 & highlevel.dataframe$T==1),]
  
  cleaned.highlevel.dataframe<-rbind(cleaned.highlevel.dataframe, highlevel.dataframe[which(highlevel.dataframe$x.counts!=5),])
  highlevel.dataframe<-highlevel.dataframe[-which(highlevel.dataframe$x.counts!=5),]
  
  for(row.index in sequence(dim(highlevel.dataframe)[1])) {
    print(paste("row.index ",row.index, " of ",dim(highlevel.dataframe)[1]))
    this.focal<-highlevel.dataframe$focal[row.index]
    this.T<-highlevel.dataframe$T[row.index]
    this.D<-highlevel.dataframe$D[row.index]
    this.x.counts<-highlevel.dataframe$x.counts[row.index]
    valid<-TRUE
    if (this.x.counts==5 & (this.D>=4 | this.D==2)) {
      
    submatrix<-highlevel.dataframe[(row.index+1):dim(highlevel.dataframe)[1],]
    submatrix<-submatrix[which(submatrix$T==this.T & submatrix$D==this.D & submatrix$x.counts==this.x.counts ),]
    for (next.row.index in sequence(dim(submatrix)[1])) {
      next.focal<-submatrix$focal[next.row.index]
        
       if(paste(grepl("x",strsplit(this.focal,"")[[1]]),collapse="")==paste(grepl("x",strsplit(next.focal,"")[[1]]),collapse="")) {
          valid<-FALSE
           print(rbind(highlevel.dataframe[row.index,1:10],submatrix[next.row.index, 1:10]))
          break()
         }
}
    }
    if(valid) {
      cleaned.highlevel.dataframe<-rbind(cleaned.highlevel.dataframe, highlevel.dataframe[row.index,])
    }
  }
return(cleaned.highlevel.dataframe)
}

cleaned.highlevel.dataframe<-CleanDuplicateModels(highlevel.dataframe)
save(cleaned.highlevel.dataframe, file="all.results.cleaned.Rsave", compress=TRUE)

highlevel.dataframe<-cleaned.highlevel.dataframe

#need to recalculate AIC weight
highlevel.dataframe$deltaAIC<-highlevel.dataframe$AIC-min(highlevel.dataframe$AIC)
expAIC<-exp(-0.5*highlevel.dataframe$deltaAIC)
highlevel.dataframe$AICweight<-expAIC/sum(expAIC)
highlevel.dataframe<-highlevel.dataframe[order(highlevel.dataframe$AICweight, decreasing=TRUE),]
save(highlevel.dataframe, file="all.results.cleaned.Rsave", compress=TRUE)


length(which(x.counts==5))
sum(highlevel.dataframe[which(x.counts==5),]$AICweight)

sum(highlevel.dataframe[which(highlevel.dataframe$D==1),]$AICweight)+sum(highlevel.dataframe[which(highlevel.dataframe$D==3),]$AICweight)

sum(highlevel.dataframe[which(highlevel.dataframe$T==1),]$AICweight)

highlevel.dataframe[order(highlevel.dataframe$AICweight, decreasing=TRUE)[1:12],1:10]

sum(highlevel.dataframe[which(highlevel.dataframe$focal=="xx1xxx"),]$AICweight)

sum(highlevel.dataframe[which(x.counts<5),]$AICweight)


potential.duplicates<-highlevel.dataframe[duplicated(highlevel.dataframe$lnL),]
