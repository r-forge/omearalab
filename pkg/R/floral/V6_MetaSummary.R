setwd("/Users/bomeara/Sites/RunsJan2012/Summaries")
load("Highlevel.dataframe.withrates.Rsave")
highlevel.dataframe<-highlevel.dataframe[,c(-3,-4)] #dropping misleading columns
best<-which.min(highlevel.dataframe$AIC)
print(highlevel.dataframe[best,1:13])
#plot(highlevel.dataframe$k_all,highlevel.dataframe$deltaAIC+1,log="y")
print(dim(highlevel.dataframe))
print("final dim should be 21126  4173")
print(highlevel.dataframe[which(highlevel.dataframe$deltaAIC<10),1:13])


print(unique(c(highlevel.dataframe[best,14:dim(highlevel.dataframe)[2] ])))

printUniqueValues<-function(focalRow) {
   combo<-focalRow$focal
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
   print(c(combo,combo.xto0, combo.mostxto0, combo.mismatch, combo.mismatch.xto0, combo.mismatch.mostxto0))
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
   print(rates)
   return(rates)
}

bestrows<-which(highlevel.dataframe$deltaAIC<10)
bestsummary.df<-data.frame(c(highlevel.dataframe[best,1:13],printUniqueValues(highlevel.dataframe[best,])),stringsAsFactors=FALSE)
for (i in sequence(length(bestrows))) {
  print(i)
 try(bestsummary.df<-rbind(bestsummary.df,(c(highlevel.dataframe[bestrows[i],1:13],printUniqueValues(highlevel.dataframe[bestrows[i],])))))
}
bestsummary.df<-bestsummary.df[-1,]
bestsummary.df
write.csv(bestsummary.df,file="bestsummary.csv")