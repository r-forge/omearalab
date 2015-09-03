library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_UtilityFns.R")
library(doMC)
library(foreach)
S=6
maxStringLength<-nchar(2^S)
registerDoMC(3) #This has a lot of I/O and memory, so make it run on fewer than the available number of processsors

setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")

scale.factor.best = 1.758664 #from /Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/full_bd_20000000_rescale_ntax.old.scale_11043
files <- rev(system(paste("ls -1 | grep Checkpoint | grep RSave | grep -v T | grep ", scale.factor.best, sep=""), intern=TRUE))


data.conversions <- gsub("x", "0", c("0x00xx", "0x01xx", "0x10xx", "0x11xx", "1x00xx", "1x01xx", "1x10xx", "1x11xx"))

focal.labels <- c("0x11xx", "xx1xxx", "1xxxxx", "xxx1xx", "0x1xxx", "0xx1xx", "xx11xx", "1x11xx", "0x01xx", "0x10xx")

T.vector <- sequence(5)
D.vector <- sequence(6)

load("~/Dropbox/SummaryRaw.RSave")
for (charStateI in 1:((2^S))) { 
	binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
	iLabelShort<-sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI)
	iLabelLong<-vectorToString(binaryStateIVector)
	names(summary.dataframe)[which(names(summary.dataframe) == paste("lambda",iLabelShort,sep="",collapse=""))]<-paste("lambda",iLabelLong,sep="",collapse="")
	names(summary.dataframe)[which(names(summary.dataframe) == paste("mu",iLabelShort,sep="",collapse=""))]<-paste("mu",iLabelLong,sep="",collapse="")
			print(paste("changing names for ",iLabelLong))
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
