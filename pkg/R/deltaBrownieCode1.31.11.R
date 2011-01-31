library(phylobase)
library(RBrownie) 


tree<-rcoal(10)  #make random ape phylogeny


#for junk data
tree1<-as(tree, "phylo4d") #turn ape phylogeny into phylo4d class
brownieTree<-as(tree1, "brownie") #turn phy into a brownie object
data = runif((nTips(brownieTree)), 1, 10) #generate junk data for the browniePhy
tdata(brownieTree)










iterateNonCensored<-function (phy, data) {
	phy<-as(phy, "phylo4d") #turn ape phylogeny into phylo4d class
	browniePhy<-as(phy, "brownie") #turn phy into a brownie object
	browniePhyNew<-addData(browniePhy, tip.data=data, dataTypes=contData())  #add data to browniePhy
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
	test.results.bulk<-c() #create an empty vector as a repository for the iteration over each taxset of the test.results of the runNonCensored function
	all.test.results<-data.frame() #create an empty data frame as a repository of the final results
	taxset.labels<-c(innerNodes, "all") #the runNonCensored() can only be fed character strings (i.e. taxset = "all"), therefore need to loop into function with these character strings
		for (j in 1:length(taxset.labels)){   #this for-loop runs the runNonCensored test over each TAXSET of the tree
			test.results.bulk<-runNonCensored(browniePhyNew, models=brownie.models.continuous()[1], taxsets = taxset.labels[j],treeloop=F, charloop=F)
			all.test.results<-rbind(all.test.results, test.results.bulk)
		}
		return(browniePhyNew)
		return(all.test.results)
}