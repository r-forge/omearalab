library(phylobase)
library(RBrownie) 

phy_tree<-rcoal(10)  #make random ape phylogeny

plot(phy_tree) #plot phylogeny

phy<-as(phy_tree, "phylo4d") #turn ape phylogeny into phylo4d class

browniePhy<-as(phy, "brownie") #turn phy into a brownie object

junk2 = runif((nTips(browniePhy)), 1, 10) #generate junk data for the browniePhy

browniePhyNew<-addData(browniePhy, tip.data=junk2, dataTypes=contData())  #add data to browniePhy

tdata(browniePhyNew) #check that the data is associated with the phylogeny

browniePhyNew@order<-"unknown"  #patch so descendants() for the root node works


###for when we want to index all clade combos###

innerNodes<-grep("internal", nodeType(browniePhyNew))  #pulls out the nodes that are internal

allDec<-vector("list", length(innerNodes)+1)  #creates list with as many places as there are innerNodes plus 1 to make space for the TAXSET_all taxon set

names(allDec)<-paste("TAXSET_", innerNodes,sep="") #all names are added the "TAXSET_" necssary for the brownie object to recognize them

	

for (ii in 1:length(allDec)){		#this for-loop creates the required (by RBrownie) TAXSET_all taxset; thank you Conrad for patching the rootNode mess!
	if (is.na(names(allDec[ii]))) {
	names(allDec)[ii]<-paste("TAXSET_all")
	allDec[[ii]]<-descendants(browniePhyNew, rootNode(browniePhyNew), type=c("tips"))
	taxasets(browniePhyNew, taxnames=names(allDec[ii]))<-labels(allDec[[ii]])
	}
}


for (i in 1:length(innerNodes)) {     #this for-loop creates the required (by RBrownie) TAXSETS_ for each inner node of the tree
	allDec[[i]]<-descendants(browniePhyNew, innerNodes[i], type=c("tips"))
	taxasets(browniePhyNew, taxnames=names(allDec[i]))<-labels(allDec[[i]])
}


taxasets(browniePhyNew) #what are the taxasets

plot.taxaset(browniePhyNew, 2) #plot a taxaset in the tree

test.results.bulk<-c() #create an empty vector as a repository for the iteration over each taxset of the test.results of the runNonCensored function

all.test.results<-data.frame() #create an empty data frame as a repository of the final results

taxset.labels<-c(innerNodes, "all") #the runNonCensored() can only be fed character strings (i.e. taxset = "all"), therefore need to loop into function with these character strings


for (j in 1:length(taxset.labels)){   #this for-loop runs the runNonCensored test over each TAXSET of the tree
	test.results.bulk<-runNonCensored(browniePhyNew, models=brownie.models.continuous()[1], taxsets = taxset.labels[j],treeloop=F, charloop=F)
	all.test.results<-rbind(all.test.results, test.results.bulk)
}

