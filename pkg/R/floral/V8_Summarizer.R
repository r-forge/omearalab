library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_UtilityFns.R")
library(doMC)
library(foreach)
S=6
registerDoMC(3) #This has a lot of I/O and memory, so make it run on fewer than the available number of processsors

setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")

scale.factor.best = 1.758664 #from /Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/full_bd_20000000_rescale_ntax.old.scale_11043
files <- rev(system(paste("ls -1 | grep Checkpoint | grep RSave | grep -v T | grep ", scale.factor.best, sep=""), intern=TRUE))


data.conversions <- gsub("x", "0", c("0x00xx", "0x01xx", "0x10xx", "0x11xx", "1x00xx", "1x01xx", "1x10xx", "1x11xx"))

focal.labels <- c("0x11xx", "xx1xxx", "1xxxxx", "xxx1xx", "0x1xxx", "0xx1xx", "xx11xx", "1x11xx", "0x01xx", "0x10xx")

T.vector <- sequence(5)
D.vector <- sequence(6)
result.df <- data.frame()

for (file.index in sequence(length(files))) {
	result.df.local <- data.frame()
	for(T.index in sequence(length(T.vector))) {
		for (D.index in sequence(length(D.vector))) {
			for(F.index in sequence(length(focal.labels))) {
				Fstring=vectorToString(convertFocalLabelToFocalVector(focal.labels[F.index], S=6, uncertainty="x"))
				file.to.load <- paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, "Condor.Rsave", sep="")
				result <- NULL
				try(load(file.to.load), silent=TRUE)
				if (!is.null(result)) {
					result.df.local <- rbind(result.df.local, data.frame(file=files[file.index], T=transitionModels$description[T.index], D=diversificationModels$description[D.index], F=focal.labels[F.index], AIC=result["AIC",1], stringsAsFactors=FALSE))
					print(tail(result.df, 1))
				}
			}
		}
	}
	result.df.local$deltaAIC <- result.df.local$AIC - min(result.df.local$AIC, na.rm=TRUE)
	result.df.local$AICweight <- exp(-0.5 * result.df.local$deltaAIC)
	result.df.local$AICweight <- result.df.local$AICweight / sum(result.df.local$AICweight)
	result.df <- rbind(result.df, result.df.local)
}

save(result.df, file="~/Dropbox/ConcatenatedResult.RSave")
write.csv(result.df, file="~/Dropbox/ConcatenatedResult.csv")


summary.dataframe <- data.frame()
for (file.index in sequence(length(files))) {
	for(T.index in sequence(length(T.vector))) {
		for (D.index in sequence(length(D.vector))) {
			print(paste(file.index, T.index, D.index))
			for(F.index in sequence(length(focal.labels))) {
				Fstring=vectorToString(convertFocalLabelToFocalVector(focal.labels[F.index], S=6, uncertainty="x"))
				file.to.load <- paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, "Condor.Rsave", sep="")
				result <- NULL
				try(load(file.to.load), silent=TRUE)
				if (!is.null(result)) {
					final.matrix.all <- result
					focalVector <- convertFocalLabelToFocalVector(focal.labels[F.index], S=6, uncertainty="x")
					qIndices<-grep("^q\\d",row.names(final.matrix.all),perl=TRUE)
					lambdaIndices<-grep("^lambda\\d",row.names(final.matrix.all),perl=TRUE)
					muIndices<-grep("^mu\\d",row.names(final.matrix.all),perl=TRUE)
					transitionModelIndex <- T.index
					diversificationModelIndex <- D.index
					tmp.dataframe<-data.frame(paste(getFocalSummaryLabel(focalVector,S=S,any="x"),sep="",collapse=""),transitionModelIndex,transitionModels[transitionModelIndex,4],diversificationModelIndex,diversificationModels[diversificationModelIndex,5],final.matrix.all[which(row.names(final.matrix.all)=="lnLik"),1],final.matrix.all[which(row.names(final.matrix.all)=="AIC"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_all"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_q"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_lambda"),1],final.matrix.all[which(row.names(final.matrix.all)=="k_mu"),1]) 
					names(tmp.dataframe)<-c("focal","T","TransitionModel","D","DiversificationModel","lnLik","AIC","k_all","k_q","k_lambda","k_mu")	
					tmp.dataframe<-cbind(tmp.dataframe,data.frame(matrix(final.matrix.all[qIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[qIndices,1])))),data.frame(matrix(final.matrix.all[lambdaIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[lambdaIndices,1])))),data.frame(matrix(final.matrix.all[muIndices,1],nrow=1,dimnames=list("",names(final.matrix.all[muIndices,1])))))
					tmp.dataframe <- cbind(tmp.dataframe, data.frame(file=files[file.index], stringsAsFactors=FALSE))
					summary.dataframe<-rbind(summary.dataframe,tmp.dataframe)
					
					
					
					
				}
			}
			save(summary.dataframe, file="~/Dropbox/SummaryRaw.RSave")
			write.csv(summary.dataframe, file="~/Dropbox/SummaryRaw.csv")

		}
	}
}

save(summary.dataframe, file="~/Dropbox/SummaryRaw.RSave")
write.csv(summary.dataframe, file="~/Dropbox/SummaryRaw.csv")
for (charStateI in 1:((2^S))) { 
	binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
	iLabelShort<-sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI)
	iLabelLong<-vectorToString(binaryStateIVector)
	names(summary.dataframe)[which(names(summary.dataframe) == paste("lambda",iLabelShort,sep="",collapse=""))]<-paste("lambda",iLabelLong,sep="",collapse="")
	names(summary.dataframe)[which(names(summary.dataframe) == paste("mu",iLabelShort,sep="",collapse=""))]<-paste("mu",iLabelLong,sep="",collapse="")
#			print(paste("changing names for ",iLabelLong))
	for (charStateJ in 1:((2^S))) { 
		binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
		jLabelShort<-sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ)
		jLabelLong<-vectorToString(binaryStateJVector)
		numberMismatches=vectorMismatch(binaryStateIVector,binaryStateJVector) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
		if (numberMismatches==1) {
			names(summary.dataframe)[which(names(summary.dataframe) == paste("q",iLabelShort,jLabelShort,sep="",collapse=""))]<-paste("q",iLabelLong,"_",jLabelLong,sep="",collapse="")			
		}
		else {
			names(summary.dataframe)[which(names(summary.dataframe) == paste("q",iLabelShort,jLabelShort,sep="",collapse=""))]<-paste("q",iLabelLong,"_",jLabelLong,"_disallowed",sep="",collapse="")			
		}
	}
}

save(summary.dataframe, file="~/Dropbox/SummaryPretty.RSave")
write.csv(summary.dataframe, file="~/Dropbox/SummaryPretty.csv")
