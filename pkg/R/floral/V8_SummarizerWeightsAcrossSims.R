load("~/Dropbox/SummaryWithWeightsAndAverages.RSave")
summary.with.weights$model.name <- paste(summary.with.weights$focal, "_T", summary.with.weights$T, "_", summary.with.weights$TransitionModel, "_D", summary.with.weights$D, "_", summary.with.weights$DiversificationModel, sep="")
all.models <- unique(summary.with.weights$model.name)
all.sims <- unique(summary.with.weights$file)
model.summaries <- data.frame(matrix(nrow=length(all.sims), ncol=length(all.models)))
for (i in sequence(length(all.models))) {
	print(i/length(all.models))
	for (j in sequence(length(all.sims))) {
		focal.data <- subset(summary.with.weights, model.name==all.models[i] & file==all.sims[j])
		if(dim(focal.data)[1]>0) {
			model.summaries[j,i] <- focal.data$AICweight	
		}	
	}
}
colnames(model.summaries) <- all.models
rownames(model.summaries) <- all.sims
save(model.summaries, file="~/Dropbox/SummaryOfWeights.RSave")
model.summaries.ordered <- model.summaries[,order(apply(model.summaries, 2, median, na.rm=TRUE), decreasing=TRUE)]
pdf(file="~/Dropbox/SummaryOfWeightsAll.pdf", width=10, height=5)
plot(x=c(1,dim(model.summaries.ordered)[2]), y=range(model.summaries.ordered, na.rm=TRUE), log="y", xlab="Model rank", ylab="Akaike weight", type="n", bty="n")
for (i in sequence(dim(model.summaries.ordered)[2])) {
	lines(x=sequence(dim(model.summaries.ordered)[2]), y=model.summaries.ordered[i,], col=rgb(0,0,0,0.3))
}
dev.off()
system("open ~/Dropbox/SummaryOfWeightsAll.pdf")
save(model.summaries.ordered, file="~/Dropbox/SummaryOfWeights.RSave")

cumulative.weights <- t(apply(model.summaries.ordered, 1, cumsum))
cumulative.weights[which(is.na(cumulative.weights))] <- 1

pdf(file="~/Dropbox/SummaryOfWeightsCumulative.pdf", width=10, height=5)
plot(x=c(1,dim(cumulative.weights)[2]), y=range(cumulative.weights, na.rm=TRUE), log="", xlab="Model rank", ylab="Cumulative Akaike weight", type="n", bty="n")
for (i in sequence(dim(cumulative.weights)[2])) {
	lines(x=sequence(dim(cumulative.weights)[2]), y= cumulative.weights[i,], col=rgb(0,0,0,0.3))
}
dev.off()
system("open ~/Dropbox/SummaryOfWeightsCumulative.pdf")

all.T <- unique(summary.with.weights$T)
model.summaries.by.T <- data.frame(matrix(nrow=length(all.sims), ncol=length(all.T)))
for (i in sequence(length(all.T))) {
	for (j in sequence(length(all.sims))) {
		focal.data <- subset(summary.with.weights, T==all.T[i] & file==all.sims[j])
		if(dim(focal.data)[1]>0) {
			model.summaries.by.T[j,i] <- sum(focal.data$AICweight, na.rm=TRUE)
		}	
	}
}
colnames(model.summaries.by.T) <- all.T
rownames(model.summaries.by.T) <- all.sims

all.D <- unique(summary.with.weights$D)
model.summaries.by.D <- data.frame(matrix(nrow=length(all.sims), ncol=length(all.D)))
for (i in sequence(length(all.D))) {
	for (j in sequence(length(all.sims))) {
		focal.data <- subset(summary.with.weights, D==all.D[i] & file==all.sims[j])
		if(dim(focal.data)[1]>0) {
			model.summaries.by.D[j,i] <- sum(focal.data$AICweight, na.rm=TRUE)
		}	
	}
}
colnames(model.summaries.by.D) <- all.D
rownames(model.summaries.by.D) <- all.sims

all.focal <- unique(summary.with.weights$focal)
model.summaries.by.focal <- data.frame(matrix(nrow=length(all.sims), ncol=length(all.focal)))
for (i in sequence(length(all.focal))) {
	for (j in sequence(length(all.sims))) {
		focal.data <- subset(summary.with.weights, focal==all.focal[i] & file==all.sims[j])
		if(dim(focal.data)[1]>0) {
			model.summaries.by.focal[j,i] <- sum(focal.data$AICweight, na.rm=TRUE)
		}	
	}
}
colnames(model.summaries.by.focal) <- all.focal
rownames(model.summaries.by.focal) <- all.sims

model.combos.weights <- array(dim=c(length(all.T), length(all.D), length(all.sims)))
for (t in sequence(length(all.T))) {
	for (d in sequence(length(all.D))) {
		for (s in sequence(length(all.sims))){
			focal.df <- subset(summary.with.weights, T==all.T[t] & D==all.D[d] & file==all.sims[s])
			model.combos.weights[t,d,s] <- sum(focal.df$AICweight, na.rm=TRUE)
		}
	}
}
			
