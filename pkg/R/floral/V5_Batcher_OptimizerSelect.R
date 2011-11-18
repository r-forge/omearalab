library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
source("/data/abc/RunsNov2011/UnifiedApproachScripts/V5_UtilityFns.R")
find.mle.method="optim"


focalVectorList<-getAllInterestingFocalVectorsStringsEfficient(S)

totalRuns<-0
runsInFile<-0
pbsCommands=""

for (focalIndex in 1:length(focalVectorList)) {
	focalVector<-stringToVector(unlist(focalVectorList[[focalIndex]]))
	for (transitionModelIndex in 1:dim(transitionModels)[1]) {
		mkdirCmd=paste("mkdir -p ",paste("/data/abc/RunsNov2011/ActualRuns",find.mle.method,"/T",transitionModelIndex,sep="",collapse=""),sep="",collapse="")
		suppressWarnings(system(mkdirCmd))
		for (diversificationModelIndex in 1:dim(diversificationModels)[1]) {
			mkdirCmd=paste("mkdir -p ",paste("/data/abc/RunsNov2011/ActualRuns",find.mle.method,"/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,sep="",collapse=""),sep="",collapse="")
			suppressWarnings(system(mkdirCmd))
			if (numberFocalCombos(focalVector) >= transitionModels$min_focalcombos[transitionModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
				if(numberFocalCombos(focalVector) >= diversificationModels$min_focalcombos[diversificationModelIndex]) { #if there aren't enough combos to make the model appropriate, don't run it
					totalRuns<-totalRuns+1
					#yay! Now we can run!
					nameRoot<-paste("T",transitionModelIndex,"_D",diversificationModelIndex,"_",vectorToString(getFocalSummaryLabel(focalVector,S,"x")),sep="",collapse="")
					dirRoot<-paste("/data/abc/RunsNov2011/ActualRuns",find.mle.method,"/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,"/",nameRoot,sep="",collapse="")
					mkdirCmd=paste("mkdir -p ",dirRoot,sep="",collapse="")
					suppressWarnings(system(mkdirCmd))
					lsString=paste(paste("ls -1 ",dirRoot,' | grep -c final.matrix.all',sep="",collapse=""))
					#print(lsString)
					finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
					if(finalMatrixAllCount==0) {
						runCommand=paste("source('/data/abc/RunsNov2011/UnifiedApproachScripts/V5_Commands_OptimizerSelect.R')\ntry(doUnifiedRun(F='",vectorToString(focalVector),"',T=",transitionModelIndex,",D=",diversificationModelIndex,",S=",partitionSize,", find.mle.method='",find.mle.method,"'))",sep="",collapse="")
						cat(runCommand,file=paste(dirRoot,'/run.R',sep=""),append=FALSE)
						if (runsInFile==0) {
							pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
							#queue="long*"
							#if (partitionSize==1) {
							#	queue="short*" #2 hr
							#}
							#else if (partitionSize<=3) {
							#	queue="medium*" #24 hr
							#}
							queue="long*"
							#if (runif(1,0,1)<0.5) { #half the time put in long queue
							#	queue="long*"
							#}
							pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
							pbsCommands=paste(pbsCommands,'#$ -M omeara.brian@gmail.com', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
							pbsCommands=paste(pbsCommands,"\n","#$ -N T",transitionModelIndex,"D",diversificationModelIndex,"F",vectorToString(getFocalSummaryLabel(focalVector,S,"x")),sep="")
						}
						pbsCommands=paste(pbsCommands,"\n",'cd /data/abc/RunsNov2011/ActualRuns',find.mle.method,'/T',transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,"/",nameRoot,sep="",collapse="")
						pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH run.R",sep="")
						#print(paste(paste("../ActualRuns",find.mle.method,"/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/run.sh',sep=""))
						pbsCommands=paste(pbsCommands,"\nrm ",' *.csv *.t ',sep="")
						runsInFile=runsInFile+1
						print(paste("Queuing run ",nameRoot," at ",date(),sep="",collapse=""))
						if (runsInFile>6) { #change this to deal with remnants
							cat(pbsCommands,file=paste(dirRoot,'/run.sh',sep=""),append=FALSE)
							#cat(pbsCommands,file=paste('/usr/bin/tail -40 ', dirRoot,'/run.Rout > ',dirRoot,'/tail.run.Rout',sep=""),append=FALSE)
							#cat(pbsCommands,file=paste('/usr/bin/zip ', dirRoot,'/run.zip  ',dirRoot,'/run.Rout',sep=""),append=FALSE)
							#cat(pbsCommands,file=paste('/usr/bin/rm ', dirRoot,'/run.Rout',sep=""),append=FALSE)
							
							print(pbsCommands)
							#print(paste("cd ",paste("../ActualRuns",find.mle.method,"/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,sep=""))
							origWD<-getwd()
							setwd(paste(paste("/data/abc/RunsNov2011/ActualRuns",find.mle.method,"/T",transitionModelIndex,"/T",transitionModelIndex,"_D",diversificationModelIndex,sep="",collapse=""),"/",nameRoot,sep=""))
							system("pwd")
							system("chmod u+x run.sh")
							system("/opt/sge/bin/lx24-amd64/qsub run.sh")
							setwd(origWD)
							Sys.sleep(10)
							runsInFile=0
							pbsCommands=""
						}
						while(as.numeric(system("/opt/sge/bin/lx24-amd64/qstat | grep -c bomeara",intern=TRUE))>500) {
							Sys.sleep(37)
						}
					}			
				}
			}
		}
	}
}