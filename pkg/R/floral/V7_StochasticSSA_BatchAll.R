setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsApril2013")
system("cp /Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7* .")
constraint.vector <- c("full", "transonly", "divonly", "symmetry")
net.div.vector <- c(TRUE, FALSE)
for (i in sequence(length(net.div.vector))) {
	for (j in sequence(length(constraint.vector))) {
		MakeRunFiles(constraint=constraint.vector[j], net.div=net.div.vector[i], submit=TRUE)
	}
}