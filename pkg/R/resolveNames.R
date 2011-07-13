library(rjson)
#takes a list of names and sends it to the iPlant TNRS site (http://tnrs.iplantcollaborative.org/)
#names<-c("zea mays","acacia","solanum","saltea","rosa_rugoso")
#returnedNames<-resolveNames(names)
#print(returnedNames)

resolveNames<-function(names,maxPerCall=Inf,verbose=TRUE) {
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
         if(verbose) {
            print(paste("finished ",nameIndex,"of ",length(names),"names")) 
         }
         startingPosition<-nameIndex+1
         namesInCall<-0
         callActual<-callBase
     }
  }
  print("Ignore a warning message about incomplete final line")
  if (length(newNames)!=length(names)) {
    warning(paste("the input name list was",length(names),"long but the new one is ",length(newNames),"long"))
  }
  return(newNames)
}