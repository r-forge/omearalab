load("~/Dropbox/SummaryRaw.RSave")
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_UtilityFns.R")
maxStringLength <- nchar(2^S)
files <- unique(summary.dataframe$file)
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

DoWeightedMean <- function(x, weights) {
	x.try <- as.numeric(x)	
	if(is.na(x.try[1])) {
		return(x[1])	
	}
	results <- NA
	try(results <- weighted.mean(x.try, w=weights))
	return(results)
}
summary.of.summary <-list()
overall.names <- c()
for (file.index in sequence(length(files))) {
	local <- subset(summary.dataframe, file==files[file.index])	
	local$deltaAIC <- local$AIC - min(local$AIC, na.rm=TRUE)
	local$AICweight <- exp(-0.5 * local$deltaAIC)
	local$AICweight <- local$AICweight / sum(local$AICweight)
	overall.names <- colnames(local)
	local.mean <- apply(local, 2, DoWeightedMean, weights=local$AICweight )
	print(local.mean)
	summary.of.summary[[file.index]]  <- local.mean
	#colnames(summary.of.summary) <- names(local.mean)
}

#colnames(summary.of.summary) <- colnames(summary.dataframe)
summary2<-data.frame(matrix(nrow=length(summary.of.summary), ncol=0))
for (i in sequence(dim(summary.dataframe)[2])) {
	values <- c()
	for (j in sequence(length(summary.of.summary))) {
		values <- append(values, summary.of.summary[[j]][i])	
	}
	if(grepl("q",colnames(summary.dataframe)[i]	) | grepl("lambda",colnames(summary.dataframe)[i]) | grepl("mu",colnames(summary.dataframe)[i]	) ) {
		values <- as.numeric(as.character(values))
		summary2 <- cbind(summary2, as.data.frame(matrix(values, nrow=1)))
	} else {
		summary2 <- cbind(summary2, as.data.frame(matrix(values, nrow=1)))
	}
}
#colnames(summary2) <- colnames(summary.dataframe)
