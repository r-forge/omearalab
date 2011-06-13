setwd("/Users/bomeara/Dropbox/InProgress/Carstens/Basic")
popVector<-rep(7,3)
migrationArray<-generateMigrationIndividuals(popVector)
createAssignment(popVector)
allValues<-c(20,.2,.02,.0001,.001,.01)
names(allValues)<-c("collapse_1","n0multiplier_1","n0multiplier_2","migration_1","migration_2","migration_3")
print(allValues)
trueModel<-57
#trueModel<-1
library(phyclust) #just to call ms
numModels<-length(migrationArray)
saveMS(popVector,migrationArray[[trueModel]],allValues,nTrees=25,file="obs.tre")
for (model in 1:numModels) {
  saveMS(popVector,migrationArray[[model]],allValues,nTrees=10000,file="sim.tre")
  output<-system("perl compareclades.pl -aassign.txt -oobs.tre -ssim.tre",intern=TRUE)
  print(paste(model,sum(as.numeric(output))))
}