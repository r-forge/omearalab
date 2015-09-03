load("RatesRaw.Rdata")
library(hisse)
library(doMC)
library(parallel)
registerDoMC(cores=detectCores())

scale.factor.best = 1.758664 #from /Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/full_bd_20000000_rescale_ntax.old.scale_11043
scale.factor.vector <- seq(from=1, to=scale.factor.best, length.out=4)


DoSingleRun <- function(x, lambda.means, mu.means, q.means, t.rescale, max.taxa, checkpoint.frequency=100, x0=1, nstart=2, make.yule=TRUE, scale.factor.vector) {
	scale.factor <- sample(scale.factor.vector,1)
	lambda.means <- scale.factor * lambda.means
	names(lambda.means) <- gsub("lambda", "", names(lambda.means))
	mu.means <- scale.factor * mu.means
	names(mu.means) <- gsub("mu", "", names(mu.means))

	if(make.yule) {
		lambda.means <- lambda.means - mu.means
		mu.means <- 0*mu.means
	}
	turnover=lambda.means + mu.means
	ef=mu.means/lambda.means
	rate.mat <- matrix(0, nrow=length(lambda.means), ncol=length(lambda.means))
	rownames(rate.mat) <- gsub("lambda", "", names(lambda.means))
	colnames(rate.mat) <- gsub("lambda", "", names(lambda.means))
	for (i in sequence(length(lambda.means))) {
	  for (j in sequence(length(lambda.means))) {
		rate.name <- paste("q",rownames(rate.mat)[i], "_", colnames(rate.mat)[j], sep="")
		matching.element <- which(names(q.means)==rate.name)
		if(length(matching.element)==1) {
			rate.mat[i,j] <- q.means[matching.element]	
		}
	  }
	}
	diag(rate.mat) <- NA

	result.sim <- NULL
	rep.id <- x
	try(result.sim <- SimulateHisse(turnover.rates=turnover, eps.values=ef, transition.rates=rate.mat, max.t=t.rescale, max.taxa=5e6, checkpoint.file=paste("Checkpoint_", rep.id, "_scale.factor_", round(scale.factor,6), "_.RSave", sep=""), checkpoint.frequency=100, x0=1, nstart=2))
	if(is.null(result.sim)) {
		return(NA)
	 } else {
		save(result.sim, file=paste("Simulation_", rep.id, "_scale.factor_", round(scale.factor,6), "_Done.RSave", sep=""))
		return(result.sim$n.surviving)
	}
}

foreach(x=sequence(100)) %dopar% DoSingleRun(x, lambda.means, mu.means, q.means, t.rescale=t.rescale, max.taxa=5e6, scale.factor.vector=scale.factor.vector)

# try(result.sim <- SimulateHisse(turnover.rates=turnover, eps.values=ef, transition.rates=rate.mat, max.t=t.rescale, max.taxa=5e6, checkpoint.file=paste("Checkpoint_", rep.id, "_.RSave", sep=""), checkpoint.frequency=100, x0=1, nstart=2))
# if(!is.null(result.sim)) {
	# save(result.sim, file=paste("Simulation_", rep.id, "_Done.RSave", sep=""))
	# num.surviving <- append(num.surviving, result.sim$n.surviving)
	# print(quantile(num.surviving))
# }
# }
