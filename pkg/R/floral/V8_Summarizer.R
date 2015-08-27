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
