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

# Definitions:
#   character: single trait, like petal symmetry
#	character combination: combination of traits, like 0100011
#   focal set: group of character combinations, like c(0100011, 0100010) is one possible focal set

# max number of parameters in this case is just 8, not bad for 500 taxa. Want to only do one case where all states are equal for all rates

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

#comboDecimal go from 1:2^S
comboAsBinaryVector<-function(comboDecimal,S) {
	return(digitsBase(comboState-1,ndigits=S)[,1])
}

comboAsBinaryString<-function(comboDecimal,S) {
	return(paste(comboAsBinaryVector(comboState,S),sep="",collapse=""))
}

comboAsDecimal<-function(comboBinary,S) {
	if (length(comboBinary)==1) { #we have a string, first make it a vector
		comboBinary<-strsplit(comboBinary,split="")[[1]]
	}
	return(1+todec(as.numeric(comboBinary)))
}

#focalAsBinaryVector goes from 000......000 (nothing is special) to 111......111 (all is special). Is of length 2^7.
#so if it is 110000...000  only combo numbers 1 and 2 are in the focal set
focalAsBinaryVector<-function(focalDecimal,S) {
	return(digitsBase(focalDecimal,ndigits=2^S)[,1])
}

focalAsBinaryString<-function(focalDecimal,S) {
	return(paste(focalAsBinaryVector(focalDecimal,S),sep="",collapse=""))
}

focalAsDecimal<-function(focalBinary,S) {
	if (length(focalBinary)==1) { #we have a string, first make it a vector
		focalBinary<-strsplit(focalBinary,split="")[[1]]
	}
	return(todec(as.numeric(focalBinary)))
}

convertFocalToCombos<-function(focalBinary) {
	return(which(focalBinary==1))
}


interestingFocal<-function(focalBinary,S) {
	interestingFocal<-FALSE
	focalCombos<-convertFocalToCombos(interestingFocal)
	if (length(focalCombos)<2) {
		interestingFocal<-TRUE #is interesting because a single focal trait or zero focal traits
	}
	else {
		comboMatrix<-matrix()
		for (i in 1:length(focalCombos)) {
			if (i==1) {
				comboMatrix<-matrix(focalAsBinaryVector(focalCombos[i]),nrow=1)
			}
			else {
				comboMatrix<-rbind(comboMatrix,focalAsBinaryVector(focalCombos[i]))
			}
		}
		varVector<-c()
		for (i in 1:length(focalCombos)) {
			varVector<-c(varVector, var(focalCombos[,i]))
		}
		numberInvariantSites<-length(which(varVector==0))
		if ((2^(S-numberInvariantSites))==dim(comboMatrix)[1]) {
			interestingFocal==TRUE #imagine just 3 chars. c(010,011) works as 01*, but c(010,111) does not (it is just *1*, but only part of *1*, omitting 110 and 011). This tests that. 
		}
	}
	return(interestingFocal)
}
