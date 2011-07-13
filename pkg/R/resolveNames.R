library(rjson)
#takes a list of names and sends it to the iPlant TNRS site (http://tnrs.iplantcollaborative.org/)
#b<-fromJSON(file="http://tnrs.iplantc.org/tnrsm-svc/matchNames?retrieve=best&names=zea%20mays,acacia,solanum,saltea")
resolveNames<-function(names,maxPerCall=Inf) {
  names<-sapply(names,sub,pattern="_",replacement=" ")
  names<-sapply(names,URLencode)
  callBase<-'http://tnrs.iplantc.org/tnrsm-svc/matchNames?retrieve=best&names='
  newNames<-rep(NA,length(names))
  namesInCall<-0
  callActual<-callBase
  startingPosition<-1
  for (nameIndex in sequence(length(names))) {
     namesInCall<-namesInCall+1
     callActual<-paste(callActual,names[nameIndex],",",sep="")
     if (namesInCall==maxPerCall || nameIndex==length(names)) {
         returnedValues<-fromJSON(file=callActual)$items
         for (returnIndex in sequence(length(returnedValues))) {
            newNames[startingPosition+returnIndex-1]<-returnedValues[[returnIndex]]$nameScientific 
         }
         startingPosition<-nameIndex+1
         namesInCall<-0
         callActual<-callBase
     }
  }
  print("Ignore a warning message about incomplete final line")
  return(newNames)
}

names<-c("zea mays","acacia","solanum","saltea","rosa_rugoso")
