wdt<-function() { #function to get the whole damn tree of life
  
  
}

downloadSetFromTreeBase<-function(study.id,time.delay=0) {
   study.id<-as.character(study.id)
   returnobject<-list(worked=FALSE,study.id=as.numeric(study.id))
   trees<-FALSE
   try(trees<-search_treebase(study.id, by="id.study"))
   if (class(trees)=="multiPhylo") {
     if(length(trees)>0) {
      meta<-FALSE
       try(meta<-metadata(study.id)[[1]])
       returnobject<-list(worked=TRUE,study.id=as.numeric(study.id), trees=trees,metadata=meta)
     }
   }
   Sys.sleep(time.delay)
   return(returnobject)
}

downloadAllTreeBase<-function(start.time=22,end.time=5,time.delay=15) {
  print("Getting list of all studies....")
  all <- search_metadata("", by="all")
 all.studies<-list()
 for (i in 1:length(all)) {
   print(paste("Doing ",i,"of ",length(all)))
   study.id<-get_study_id(all[i])
   try(all.studies<-c(all.studies,list(downloadSetFromTreeBase(study.id,time.delay))))
   save(all.studies,file="treebase.rda",compress=TRUE)
    # while(!(as.numeric(format(Sys.time(), "%H"))>=start.time || as.numeric(format(Sys.time(), "%H"))<end.time)) {
    #    Sys.sleep(60) 
    # }
 }
 return(all.studies) 
}