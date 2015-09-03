library(hisse)
library(geiger)
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")
source("V8_Checkrun.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")

scale.factor.best = 1.758664 #from /Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/full_bd_20000000_rescale_ntax.old.scale_11043
files <- system(paste("ls -1 | grep Checkpoint | grep RSave | grep ", scale.factor.best, sep=""), intern=TRUE)


data.conversions <- gsub("x", "0", c("0x00xx", "0x01xx", "0x10xx", "0x11xx", "1x00xx", "1x01xx", "1x10xx", "1x11xx"))

focal.labels <- c("0x11xx", "xx1xxx", "1xxxxx", "xxx1xx", "0x1xxx", "0xx1xx", "xx11xx", "1x11xx", "0x01xx", "0x10xx")

T.vector <- sequence(5)
D.vector <- sequence(6)
#focal.vectors <- sapply(focal.labels, convertFocalLabelToFocalVector, S=6, uncertainty="x")

for (file.index in sequence(length(files))) {
	load(files[file.index])
	phy <- SimToPhylo(checkpoint.result$results, include.extinct=FALSE, drop.stem=TRUE)
	subsetted <- treedata(drop.random(phy, Ntip(phy)-464), data=phy$tip.state, warnings=FALSE, sort=TRUE)
	phy <- subsetted$phy
	data <- subsetted$data
	data.new <- data.conversions[data+1]
	names(data.new) <- rownames(data)
	phy$tip.state <- data.new
	for(T.index in sequence(length(T.vector))) {
		for (D.index in sequence(length(D.vector))) {
			for(F.index in sequence(length(focal.labels))) {
				Fstring=vectorToString(convertFocalLabelToFocalVector(focal.labels[F.index], S=6, uncertainty="x"))
				result<-NULL
				try(result<-doSingleModelFromHisseSim(phy, F= Fstring, T=T.vector[T.index], D=D.vector[D.index], S=6))
				save(result, file.index, files, phy, T.index, D.index, Fstring, file=paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, ".Rsave", sep="")) 
			}
		}
	}
	
}
