library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions) #for converting from binary back to decimal


#basic idea:
# 7 binary traits
# means 2^7 = 128 possible character-state combinations (note: not all of them are realized)
# select between 0 and 64 of these combos as "focal" ones
# syntax for focal set: a 128-long binary string, telling whether a given combo is in the focal set
# filter these for meaningful sets: a set consisting of combos that are 0*1**** is meaningful, for example
# [problem: equal sets: 0****** vs 1****** are the same under full models (i.e., different birth, death, and transition rates -- but only in that case]

# transition models:
#	There are four transition rates q_focal->nonfocal, q_focal->focal, q_nonfocal->nonfocal, q_nonfocal->focal = qFN, qFF, qNF, qNF. Note that you can have within set transitions where the number of combos in a set is greater than 1: set 0*1**** has a 0010000->0010001 transition rate
#	1: all rates the same: qNF=qNN=qFF=qFN
#	2: inflow different: qNF vs qNN=qFF=qFN
#	3: outflow different: qFN vs qNN=qFF=qNF
#	4: inflow and outflow different: qFN vs qNF vs qFF=qNN
#	5: freedom: qFN vs qNF vs qFF vs qNN

# diversification models:
#	There are four diversification rates: bF, bN, dF, dN
#	1: yule: bF=bN, dF=dN=0
#	2: two-rate yule: bF vs bN vs dF=dN=0
#	3: simple birth-death: bF=bN, dF=dN
#	4: two birth, one death: bF vs bN vs dF=dN
#	5: one birth, two death: bF=bN vs dF vs dN
#	6: freedome: bF vs bN vs dF vs dN

# max number of parameters in this case is just 8, not bad for 500 taxa

# constrain root state = 0000000

#Global definitions
nchar=7
S=nchar

#utility functions
vectorMismatch<-function(vector1, vector2) {
	if (length(vector1)!=length(vector2)) {
		return(NA)
	}
	else {
		return(sum(1-(binaryStateIVector==binaryStateJVector))) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
	}
}

#char states go from 1:2^S
charAsBinaryVector<-function(charState,S) {
	return(digitsBase(charState-1,ndigits=S)[,1])
}

charAsBinaryString<-function(charState,S) {
	return(paste(charAsBinaryVector(charState,S),sep="",collapse=""))
}

charAsState<-function(charBinary,S) {
	if (length(charBinary)==1) { #we have a string, first make it a vector
		charBinary<-strsplit(charBinary,split="")[[1]]
	}
	return(1+todec(as.numeric(charBinary)))
}

