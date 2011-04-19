library(diversitree) #obvious
library(sfsmisc) #for counting in binary
nchar=7

numberOfTransitionModelsForPartitionSize<-c(4,17,5+3*(2^3),5+3*(2^4),5+3*(2^5),5+3*(2^6),5+3*(2^7)) #make sure to include models that reduce to independent models: though they are present in smaller model sets, we might want independent transitions but dependent diversification
numberOfDiversificationModelsForPartitionSize<-c(19,12,6+3*(2^3),6+3*(2^4),6+3*(2^5),6+3*(2^6),6+3*(2^7))


#three steps:
#Step 1
#	Get assignment of components: i.e., bisse for trait 1, musse for traits 1+2, musse for traits 3+5+6
#	For each component, try each of the possible models: i.e., for bisse, full, equal, uncorrelated; for 3 char models, maybe the "one is odd" models, split models, etc.
#   Record for each of these the model partition, combination type, model type (transition + diversification type), K, lnL, name of file containing the values
#
#Step 2
#	Get AIC of all the model combinations: i.e., AIC of aabbbcc where a=musse uniform, b=musse uncorrelated, etc.
#	Find the best model combinations
#	Report weights for different model types: how important is correlation between chars 2 and 3?
#
#Step 3
#	Get parameter estimates from best? equally good? average? models

#Step 1
#we will store these results in a dataframe
#char 1, 2, 3.... \t combosize\t model \t lnL \t K \t AIC \toutput file name
#but first have to make them
numberOfModels=0;
runsInFile=0;
pbsCommands=""
partitionIterator=(2^nchar)-1 #doing 1234567 ONLY
#print(digitsBase(partitionIterator,ndigits=nchar))
partitionScheme<-digitsBase(partitionIterator,ndigits=nchar)
partitionScheme<-partitionScheme[,1] #1 means include that char in partition, NA means do not
partitionScheme<-partitionScheme*c(1:nchar) #so this will result in a vector with digits showing if the character is included
partitionSize<-length(which(partitionScheme>0))
partitionSchemeText=paste(partitionScheme,sep="",collapse="_")
mkdirCmd=paste("mkdir ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),sep="",collapse="")
suppressWarnings(system(mkdirCmd))
for (transitionIndex in 1:numberOfTransitionModelsForPartitionSize[partitionSize]) { 
	for (diversificationIndex in 1:numberOfDiversificationModelsForPartitionSize[partitionSize]) {
		if (min(numberOfTransitionModelsForPartitionSize[partitionSize],numberOfDiversificationModelsForPartitionSize[partitionSize])>0) { #just to make sure we have models: for (1:1) and for (1:0) still run
			#Sys.sleep(1) #because I'm a nice guy
			numberOfModels=numberOfModels+1
			#print(numberOfModels)
			nameRoot=paste("P",partitionSchemeText,"_T",transitionIndex,"_D",diversificationIndex,sep="",collapse="")
			mkdirCmd=paste("mkdir ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,sep="",collapse="")
			mkdirResult<-suppressWarnings(system(mkdirCmd))
			if (mkdirResult!=0) { #already exists. We might want to rerun stuff if earlier runs failed, but we'll remove pointless files nonetheless
				system(paste("rm ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/*t',sep="",collapse=""))
				system(paste("rm ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/*csv',sep="",collapse=""))
				system(paste("rm ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/*sh',sep="",collapse=""))
				system(paste("rm ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/*sh',sep="",collapse=""))
				system(paste("rm ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/*fit.final',sep="",collapse=""))
				system(paste("rm ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/*final.matrix',sep="",collapse=""))
			}
			lsString=paste(paste("ls -1 ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,' | grep -c final.matrix.all',sep="",collapse=""))
			print(lsString)
			finalMatrixAllCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
			print(paste("finalMatrixAllCount = ",finalMatrixAllCount," for ",nameRoot))
			if(finalMatrixAllCount==0) {
				runCommand=paste("source('../../../UnifiedApproachScripts/UnifiedApproachV3Commands1234567Only.R')\ndoUnifiedRun(P='",partitionSchemeText,"',T=",transitionIndex,",D=",diversificationIndex,",S=",partitionSize,")",sep="",collapse="")
				cat(runCommand,file=paste(paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/run.R',sep=""),append=FALSE)
				if (runsInFile==0) {
					pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
					#queue="long*"
					#if (partitionSize==1) {
					#	queue="short*" #2 hr
					#}
					#else if (partitionSize<=3) {
					#	queue="medium*" #24 hr
					#}
					queue="medium*"
					pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
					pbsCommands=paste(pbsCommands,'#$ -M omeara.brian@gmail.com', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
					pbsCommands=paste(pbsCommands,"\n","#$ -N R_",gsub("_","",partitionSchemeText),"\n", 'module load R/2.12.0',sep="")
				}
				pbsCommands=paste(pbsCommands,"\n",'cd /data/abc/RunsFebruary2011/ActualRuns/P',partitionSchemeText,"/",nameRoot,sep="",collapse="")
				pbsCommands=paste(pbsCommands,"\n","/data/apps/R/R-2.12.0/bin/R CMD BATCH run.R",sep="")
				#print(paste(paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/run.sh',sep=""))
				pbsCommands=paste(pbsCommands,"\nrm ",' *.csv *.t ',sep="")
				runsInFile=runsInFile+1
				if (runsInFile>200) { #change this to deal with remnants
					cat(pbsCommands,file=paste(paste("/data/abc/RunsFebruary2011/ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,'/run.sh',sep=""),append=FALSE)
					print(pbsCommands)
					#print(paste("cd ",paste("../ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,sep=""))
					setwd(paste(paste("/data/abc/RunsFebruary2011/ActualRuns/P",partitionSchemeText,sep="",collapse=""),"/",nameRoot,sep=""))
					system("pwd")
					system("chmod u+x run.sh")
					system("qsub run.sh")
					setwd("../..")
					Sys.sleep(1)
					runsInFile=0
					pbsCommands=""
				}
			}
			while(as.numeric(system("qstat | grep -c bomeara",intern=TRUE))>1800) {
				Sys.sleep(117)
			}
		}
	}
}



