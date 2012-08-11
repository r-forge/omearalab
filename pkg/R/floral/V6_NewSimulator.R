library(RColorBrewer)

OMearaSSA<-function(x0, q.vector, lambda.vector, mu.vector, tf, maxWallTime, verbose=TRUE) {
  #make sure order of states is the same in all input objects
  q.matrix<-matrix(0,nrow=length(x0),ncol=length(x0))
  rownames(q.matrix)<-names(x0)
  colnames(q.matrix)<-names(x0)
  for (i in sequence(length(q.vector))) {
    q.matrix[ which(rownames(q.matrix)==extractFrom(names(q.vector)[i])), 
              which(colnames(q.matrix)==extractTo(names(q.vector)[i]))
              ]<-q.vector[i]
  }
  names(lambda.vector)<-gsub("lambda","",names(lambda.vector))
  names(mu.vector)<-gsub("mu","",names(mu.vector))
  #print(q.matrix)
  all.matrix<-cbind(q.matrix,matrix(lambda.vector,ncol=1),matrix(mu.vector,ncol=1))
  #print(all.matrix)
  start.time<-Sys.time()
  x<-x0
  history<-matrix(c(0,x0),nrow=1)
  while ((Sys.time()-start.time)<maxWallTime) {
    rate.matrix<-ScaleTransitionProbs(all.matrix,x)
    time.interval<-rexp(1, sum(rate.matrix))
    if ((sum(history[,1])+time.interval)>tf) {
      history<-rbind(history,c(tf-sum(history[,1]),x)) #done running
      return(history)
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

    history<-rbind(history,c(time.interval,x))
    if(verbose) {
     print(paste(sum(x),paste(x,sep=" ",collapse=" "),time.interval, sum(history[,1]), move, sep=" "))
    }
    if(sum(x)==0) {
      if(verbose) {
        print("all extinct")
      }
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
 plot(range(cumsum(history[,1])),range(history[,2:dim(history)[2] ]),xlab="time",ylab="ntaxa",type="n",bty="n") 
 for (i in 2:dim(history)[2]) {
  lines( cumsum(history[,1]), history[, i],col=mypalette[i-1])
 }
  legend(x="topleft",legend=colnames(history)[2:(ncombos+1)], text.col=mypalette)
}
  

doParallelSSA<-function(x0, q.means, lambda.means, mu.means, tf=130, maxWallTime=Inf, verbose=F) {
  file.name<-paste("SSA_",format(Sys.time(), "%b%d_%H_%M_%S"),"_",round(runif(1,1,1000000)),".RSave",sep="")
  survivors<-0
  history<-0
  attempts<-0
  while(survivors<=0) {
    attempts<-attempts+1
    history<-OMearaSSA(x0, q.means, lambda.means, mu.means, tf=tf, maxWallTime=maxWallTime, verbose=verbose)
    survivors<-history[dim(history)[1],2]
    print(c(attempts,max(apply(history[,2:9],1,sum))))
  }
  save(history, file=file.name)
}

#plot.OMearaSSA(history)
