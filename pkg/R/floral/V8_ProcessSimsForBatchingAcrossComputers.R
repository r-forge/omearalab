library(hisse)
library(geiger)
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")
source("V8_Checkrun.R")
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")

scale.factor.best = 1.758664 #from /Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/full_bd_20000000_rescale_ntax.old.scale_11043
files <- rev(system(paste("ls -1 | grep Checkpoint | grep RSave | grep -v T | grep ", scale.factor.best, sep=""), intern=TRUE))


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
				batch.file <- paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, "Batcher.Rsave", sep="")
				save(list=ls(), file=batch.file)
				batcher.file <- paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, "Batcher.R", sep="")
				prog.file <- paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, ".program", sep="")
				submit.file <- paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, ".submit", sep="")
				pkgs.needed <- c("gmp","corpcor","partitions","sfsmisc","diversitree","Rcpp","subplex","ape","deSolve" )

				cat(paste("load('",batch.file, "')\n", sep=""), file=batcher.file, append=FALSE)
				cat("\n", file=batcher.file, append=TRUE)
				cat('
					for(package.index in sequence(length(pkgs.needed))) {
						try(install.packages(pkgs.needed[package.index], repos="http://cran.us.r-project.org", lib="."))
						try(library(pkgs.needed[package.index], character.only=TRUE))
					}
				', file=batcher.file, append=TRUE)
				cat("\n", file=batcher.file, append=TRUE)
				
				cat("\ntry(result<-doSingleModelFromHisseSim(phy, F= Fstring, T=T.vector[T.index], D=D.vector[D.index], S=6))\n", file=batcher.file, append=TRUE)
				cat('save(result, file.index, files, phy, T.index, D.index, Fstring, file=paste(files[file.index],"_T",T.index, "_D",D.index, "_F_", Fstring, "Condor.Rsave", sep=""))', file=batcher.file, append=TRUE)
				cat(paste('#! /bin/sh


R CMD BATCH ',batcher.file,' 
echo "I am process id $$ on" `hostname`
', sep=""), file=prog.file)
				system(paste("chmod u+x",prog.file))
				cat(paste('executable=', prog.file,'
universe=vanilla
arguments=Example.$(Cluster).$(Process) 5
output=results.output.$(Process)
error=results.error.$(Process)
transfer_input_files=',batcher.file,',', batch.file,'
log=results.log.$(Process)
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT
queue 1

', sep=""), file=submit.file)
			system(paste("/condor/condor-installed/bin//condor_submit", submit.file))
			Sys.sleep(3)
			}
		}
	}
	
}
