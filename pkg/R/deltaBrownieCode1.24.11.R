library(phylobase)
library(RBrownie) 



phy_tree<-rcoal(10)

plot(phy_tree)

phy<-as(phy_tree, "phylo4d")

#turn phy into a brownie object

browniePhy<-as(phy, "brownie")


#generate junk data for the browniePhy

junk2 = runif((nTips(browniePhy)), 1, 10)


#add data to browniePhy

browniePhyNew<-addData(browniePhy, tip.data=junk2, dataTypes=contData())



tALL<-grep("t*", tipLabels(browniePhyNew), value=TRUE)

taxasets(browniePhyNew, taxnames="all")<-tALL


t674<-grep("t[6,7,4]", tipLabels(browniePhyNew), value=TRUE)

taxasets(browniePhyNew, taxnames="t674")<-t674


#what are the taxasets
taxasets(browniePhyNew)



#plot a taxaset in the tree
plot.taxaset(browniePhyNew, 1)

#run Non-Censored test

test.results = runNonCensored(browniePhyNew,models=brownie.models.continuous()[1:2], treeloop=F, charloop=F)
summaryCont(test.results)


#test.results = runCensored(browniePhyNew,models=brownie.models.continuous()[1:2])
test.results = runCensored(browniePhyNew, taxset = "t4t6",reps = 1000, charloop = TRUE)
summaryRatetest(test.results)


###for when we want to index all clade combos

innerNodes<-grep("internal", nodeType(browniePhyNew))

allDec<-vector("list", length(innerNodes))

names(allDec)<-paste("TAXSET_", innerNodes, sep="")

 for (i in 1:length(innerNodes)) {
	allDec[[i]]<-descendants(browniePhyNew, innerNodes[i], type=c("tips"))
 	taxasets(browniePhyNew, taxnames=names(allDec[i]))<-labels(allDec[[i]])



}
