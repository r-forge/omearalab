library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
library(doMC)
library(foreach)
source("V6_UtilityFns.R")
source("V6_Commands.R")


focalVectorList<-rev(getAllInterestingFocalVectorsStringsEfficient(S))

totalRuns<-0
runsInFile<-0
pbsCommands=""
registerDoMC(24)
for (transitionModelIndex in 1:dim(transitionModels)[1]) {
	mkdirCmd=paste("mkdir -p ",paste("../ActualRunsMac/T",transitionModelIndex,sep="",collapse=""),sep="",collapse="")
	suppressWarnings(system(mkdirCmd))
	for (diversificationModelIndex in 1:dim(diversificationModels)[1]) {
		mkdirCmd=paste("mkdir -p ",paste("../ActualRunsMac/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,sep="",collapse=""),sep="",collapse="")
		suppressWarnings(system(mkdirCmd))
	}
}

individualRun<-function(focalIndex,transitionModelIndex,diversificationModelIndex) {
	focalVector<-stringToVector(unlist(focalVectorList[[focalIndex]]))
	if (numberFocalCombos(focalVector) >= transitionModels$min_focalcombos[transitionModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
		if(numberFocalCombos(focalVector) >= diversificationModels$min_focalcombos[diversificationModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
			totalRuns<-totalRuns+1
			#yay! Now we can run!
			nameRoot<-paste("T",transitionModelIndex,"_D",diversificationModelIndex,"_",vectorToString(getFocalSummaryLabel(focalVector,S,"x")),sep="",collapse="")
			dirRoot<-paste("../ActualRunsMac/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,"/",nameRoot,sep="",collapse="")
			mkdirCmd=paste("mkdir -p ",dirRoot,sep="",collapse="")
			suppressWarnings(system(mkdirCmd))
			lsString=paste(paste("ls -1 ",dirRoot,' | grep -c final.matrix.all',sep="",collapse=""))
			#print(lsString)
			finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
			if(finalMatrixAllCount==0) {
				try(doUnifiedRun(F=vectorToString(focalVector),T=transitionModelIndex,",D=",diversificationModelIndex,",S=",partitionSize,"))",sep="",collapse="")
			}			
		}
	}
}

foreach(focalIndex=c(1:length(focalVectorList))) %:% foreach(transitionModelIndex=c(1:dim(transitionModels)[1])) %:% foreach(diversificationModelIndex=c(1:dim(diversificationModels)[1])) %dopar% { individualRun(focalIndex,transitionModelIndex,diversificationModelIndex) }