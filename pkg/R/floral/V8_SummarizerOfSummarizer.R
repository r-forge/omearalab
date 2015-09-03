library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions)
library(gmp) #for dealing with big integers
setwd("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral")
source("V6_UtilityFns.R")
S=6
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015")
scale.factor.best = 1.758664 #from /Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/full_bd_20000000_rescale_ntax.old.scale_11043


data.conversions <- gsub("x", "0", c("0x00xx", "0x01xx", "0x10xx", "0x11xx", "1x00xx", "1x01xx", "1x10xx", "1x11xx"))

focal.labels <- c("0x11xx", "xx1xxx", "1xxxxx", "xxx1xx", "0x1xxx", "0xx1xx", "xx11xx", "1x11xx", "0x01xx", "0x10xx")

T.vector <- sequence(5)
D.vector <- sequence(6)
load("~/Dropbox/SummaryPretty.RSave")
files <- unique(summary.dataframe$file)
summary.with.weights <- data.frame()
model.averages <- data.frame()
for (file.index in sequence(length(files))) {
	result.df.local <- subset(summary.dataframe, file==files[file.index])
	result.df.local$deltaAIC <- result.df.local$AIC - min(result.df.local$AIC, na.rm=TRUE)
	result.df.local$AICweight <- exp(-0.5 * result.df.local$deltaAIC)
	result.df.local$AICweight <- result.df.local$AICweight / sum(result.df.local$AICweight)
	summary.with.weights <- rbind(summary.with.weights, result.df.local)
	result.df.local.params.only <- result.df.local[,12: 4171] 
	
	model.average <- rep(NA, dim(result.df.local.params.only)[2])
	for (param.index in sequence(dim(result.df.local.params.only)[2])) {
		model.average[param.index] <- weighted.mean(result.df.local.params.only[,param.index], result.df.local$AICweight)	
	}
	model.average.df <- data.frame(matrix(model.average, nrow=1))
	colnames(model.average.df) <- colnames(result.df.local.params.only)
	model.average.df$file = files[file.index]
	model.averages <- rbind(model.averages, model.average.df)
}
save(summary.with.weights, model.averages, file="~/Dropbox/SummaryWithWeightsAndAverages.RSave")


load("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/PerformanceCheck_May2015/RatesRaw.Rdata")
lambda.means <- scale.factor.best * lambda.means
names(lambda.means) <- gsub("lambda", "", names(lambda.means))
mu.means <- scale.factor.best * mu.means
names(mu.means) <- gsub("mu", "", names(mu.means))
lambda.means <- lambda.means - mu.means
mu.means <- 0*mu.means
model.averages.disallowed.purged <- model.averages[,!grepl("disallowed", colnames(model.averages))]
div.means <- lambda.means - mu.means
names(div.means) <- paste("div", names(lambda.means), sep="")
#param.names <- gsub("x", '.', c(paste("lambda", names(lambda.means), sep=""), paste("mu", names(mu.means), sep=""), names(q.means)))
#param.true <- c(lambda.means, mu.means, q.means)
param.names <- gsub("x", '.', c(names(q.means)))
param.true <- c(q.means)
estimates.range <- range(q.means)
for(i in sequence(length(param.names))) {
	print(param.names[i]) 
	print(param.true[i])
	local.params <- model.averages.disallowed.purged[, grepl(param.names[i],colnames(model.averages.disallowed.purged))]	
	estimates <- apply(local.params, 1, mean)
	estimates.range <- range(c(estimates.range, estimates))
	print(quantile(estimates, seq(from=0, to=1, length.out=11)))
}

pdf(file="~/Dropbox/Qmeans.pdf", width=7, height=7)
plot(x=c(0.5, 0.5+length(q.means)), y=estimates.range, type="n", xlab="", ylab="transition rate", xaxt="n", bty="n", log="y")
#text(sequence(length(q.means)), rep(min(estimates.range), length(q.means)), names(q.means), cex=0.5, srt=90)
axis(side=1, at=sequence(length(q.means)), labels=names(q.means), las=2, cex.axis=0.5)
for(i in sequence(length(param.names))) {
	print(param.names[i]) 
	print(param.true[i])
	local.params <- model.averages.disallowed.purged[, grepl(param.names[i],colnames(model.averages.disallowed.purged))]	
	estimates <- apply(local.params, 1, mean)
	points(rep(i, length(estimates)), estimates, pch=20, col=rgb(0, 0, 0, 0.3))
	points(i, median(estimates), pch="-", col="purple", cex=3)
	points(i, q.means[i], pch="-", col="red", cex=2)

}
dev.off()
system("open ~/Dropbox/Qmeans.pdf")

lambda.names <- gsub("x", '.', paste("lambda", names(lambda.means), sep=""))
mu.names <- gsub("x", '.', paste("mu", names(mu.means), sep=""))
pdf(file="~/Dropbox/Div.pdf", width=7, height=7)

plot(x=c(0.5, 0.5+length(div.means)), y=c(0.04, 0.13), type="n", xlab="", ylab="diversification rate", xaxt="n", bty="n")
text(sequence(length(div.means)), rep(0.04, length(div.means)), gsub("div","",names(div.means)), pos=1, cex=0.8)
for(i in sequence(length(div.means))) {
	print(names(div.means)[i])
	print(div.means[i])
	local.params.lambda <- model.averages.disallowed.purged[, grepl(lambda.names[i],colnames(model.averages.disallowed.purged))]	
	local.params.mu <- model.averages.disallowed.purged[, grepl(mu.names[i],colnames(model.averages.disallowed.purged))]	
	local.params.div <- local.params.lambda - local.params.mu
	estimates <- apply(local.params.div, 1, mean)
	print(quantile(estimates, seq(from=0, to=1, length.out=11)))
	points(rep(i, length(estimates)), estimates, pch=20, col=rgb(0, 0, 0, 0.3))
	points(i, median(estimates), pch="-", col="purple", cex=4)
	points(i, div.means[i], pch="-", col="red", cex=3)
}

dev.off()
system("open ~/Dropbox/Div.pdf")

