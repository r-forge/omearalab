setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsApril2013")
system("cp /Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7*R .")
source("V7_StochasticSSASims_CreateAllFilesForRun.R")
constraint.vector <- c("full", "transonly", "divonly", "symmetry")
net.div.vector <- c(TRUE, FALSE)
for (i in sequence(length(net.div.vector))) {
	for (j in sequence(length(constraint.vector))) {
		MakeRunFiles(constraint=constraint.vector[j], net.div=net.div.vector[i], submit=TRUE)
	}
}


constraint.vector <- c("full")
net.div.vector <- c(TRUE, FALSE)
for (i in sequence(length(net.div.vector))) {
  for (j in sequence(length(constraint.vector))) {
    MakeRunFiles(constraint=constraint.vector[j], net.div=net.div.vector[i], x0=c(0, 2, 0, 0, 0, 0, 0, 0), submit=TRUE)
  }
}