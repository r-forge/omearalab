library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
source("/data/abc/RunsApril2011/UnifiedApproachScripts/V4_UtilityFns.R")


focalVectorList<-getAllInterestingFocalVectorsStringsEfficient(S)

totalRuns<-0
totalPossibleRuns<-0
doneRuns<-0
insufficientNumbers<-0
runsInFile<-0
pbsCommands=""
print(c("Theoretical # runs","Runs without enough combos","Possibly valid runs","Completed valid runs"))

for (focalIndex in 1:length(focalVectorList)) {
	focalVector<-stringToVector(unlist(focalVectorList[[focalIndex]]))
	for (transitionModelIndex in 1:dim(transitionModels)[1]) {
		for (diversificationModelIndex in 1:dim(diversificationModels)[1]) {
			totalPossibleRuns<-totalPossibleRuns+1
			insufficientNumbers<-insufficientNumbers+1
			if (numberFocalCombos(focalVector) >= transitionModels$min_focalcombos[transitionModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
				if(numberFocalCombos(focalVector) >= diversificationModels$min_focalcombos[diversificationModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
					insufficientNumbers<-insufficientNumbers-1 #we did run it
					totalRuns<-totalRuns+1
					#yay! Now we can run!
					nameRoot<-paste("T",transitionModelIndex,"_D",diversificationModelIndex,"_",vectorToString(getFocalSummaryLabel(focalVector,S,"x")),sep="",collapse="")
					dirRoot<-paste("/data/abc/RunsApril2011/ActualRuns/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,"/",nameRoot,sep="",collapse="")
					lsString=paste(paste("ls -1 ",dirRoot,' | grep -c final.matrix.all',sep="",collapse=""))
					finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
					lsString=paste(paste("ls -1 ",dirRoot,' | grep -c run.Rout',sep="",collapse=""))
					runRoutCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))				
					if(finalMatrixAllCount!=0) {
						if (runRoutCount==1) {
							system(paste("tail -40 ",dirRoot,"/run.Rout > ",dirRoot,"/tail.run.Rout",sep=""))
							system(paste("zip ",dirRoot,"/run.zip ",dirRoot,"/run.Rout",sep=""))
							system(paste("rm ",dirRoot,"/run.Rout",sep=""))
						}
					}
				}
			}
		}
	}
}
