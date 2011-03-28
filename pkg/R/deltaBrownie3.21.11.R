#this code should find the best break among alternate morphological models

library(phylobase)
library(RBrownie)
library(geiger) 

source("/Users/halamillo/Desktop/BrownieCode/generate.simmap.R")

#make random small ape phylogeny and generate random data for it
tree<-rcoal(10) 
phy1<-tree
data1 = runif((nTips(phy1)), 1, 10)
data1


#for three-taxa with real names example (short example)
#tip.label<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp")
#tree1<-rcoal(5, tip.label)
#tip.label2<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi")
#tree2<-rcoal(5, tip.label)


tip.label<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp")
completeTree<-rcoal(5, tip.label)

completeBipesData<-read.table("/Users/halamillo/Desktop/completeBipes.txt", row.names=1)

tip.label2<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi")
partialTree<-rcoal(4, tip.label2)

partialBipesData<-read.table("/Users/halamillo/Desktop/partialBipes.txt", row.names=1)

tip.label3<-c("Bipes_biporus", "Bipes_cannaliculatus", "Bipes_tridactylus", "Bipes_alvarezi", "Bipes_sp", "Bipes_sp2")
extraTree<-rcoal(6, tip.label3)

extraBipesData<-read.table("/Users/halamillo/Desktop/extraBipes.txt", row.names=1)


#data = runif((nTips(tree1)), 1, 10)
#data
#phy<-tree1
#bipes.results<-iterateNonCensored(tree1, data)

#for example with dipsadine tree and junk data
#phy<-read.tree("/Users/halamillo/Desktop/NonAcrochMonoRegCons2.phy")
#data<-read.table("/Users/halamillo/Desktop/dipsmorphdatatabdel.txt", row.names=1)




iterateNonCensored<-function (phy, data, name.check=TRUE) {
	if(name.check){
		
		#check that taxa match between data and tree
		checked.object<-name.check(phy,data)
	
		if (checked.object[1] == "OK") { #pillaged kindly from Banbury
			cat("Phylogeny and character matrix are in agreement...moving-on with the rate analysis!")
				
		
	
			}else {
					warning("Phylogeny and character matrix are mismatched:\n")
				cat("missmatched taxa have been dropped from analysis\n")
				#cat("(In case you are curious check your working directoy for files(s) with data that didn't match)\n")
	
				if (length(checked.object$Tree.not.data) > 0) {  
					#a<-checked.object$Tree.not.data
					#write(a, file="Tree.not.data.txt")
					phy<-drop.tip(phy, checked.object$Tree.not.data)
				}
	
				if (length(checked.object$Data.not.tree) > 0) {
					#b<-checked.object$Data.not.tree
					#write(b, file="Data.not.tree.txt")
					which(rownames(data) %in% checked.object$Data.not.tree)->rows
					as.data.frame(data[-rows,])->data1
					colnames(data1)<-colnames(data)
					rownames(data1)<-rownames(data)[-rows]
					data<-data1
				}
			}
	phy<-as(phy, "phylo4d") #turn ape phylogeny into phylo4d class
	browniePhy<-as(phy, "brownie") #turn phy into a brownie object
	brownie.tree1<-addData(browniePhy, tip.data=data, dataTypes=contData())  #add data to browniePhy
	brownie.tree1@order<-"unknown"  #patch so descendants() for the root node works
	###for when we want to index all clade combos###
	Nodes1<-grep("tip",nodeType(brownie.tree1))  #pulls out the nodes that are internal
	Nodes2<-grep("internal",nodeType(brownie.tree1))  #pulls out the nodes that are internal
	Nodes<-c(Nodes1, Nodes2)
	allDec<-vector("list", length(Nodes)+1)  #creates list with as many places as there are Nodes plus 1 to make space for the TAXSET_all taxon set
	names(allDec)<-paste("TAXSET_", Nodes,sep="") #all names are added the "TAXSET_" necessary for the brownie object to recognize them
		
		for (ii in 1:length(allDec)){		#this for-loop creates the required (by RBrownie) TAXSET_all taxset; thank you Conrad for patching the rootNode mess!
			if (is.na(names(allDec[ii]))) {
				names(allDec)[ii]<-paste("TAXSET_all")
				allDec[[ii]]<-descendants(brownie.tree1, rootNode(brownie.tree1), type=c("tips"))
				taxasets(brownie.tree1, taxnames=names(allDec[ii]))<-labels(allDec[[ii]])
			}
		}
		
		for (i in 1:length(Nodes)) {     #this for-loop creates the required (by RBrownie) TAXSETS_ for each inner node of the tree
			allDec[[i]]<-descendants(brownie.tree1, Nodes[i], type=c("tips"))
			taxasets(brownie.tree1, taxnames=names(allDec[i]))<-labels(allDec[[i]])
		}
		
	nombres<-vector("list", length(allDec))
	
		for (i in 1:length(allDec)){  #pulls out taxa that will be used in the manufacutre of taxa.vectors
			nombres[[i]]<-labels(allDec[i][[1]])
		}
		for (i in 1:length(nombres)){ # create the empty list of taxa.vectors
			names(nombres)[i]<-paste("taxa.vector.",i, sep="")
		}
	
	####so now here use each of those taxa.vectors to generate the other simmap phylogenies
	#### and read in each simmap formatted tree into an object

	
	trees.list<-c()
	trees.list.phy<-vector("list", length(nombres))
		for (i in 1:length(nombres)){
			trees.list[i]<-generate.simmap(phy, nombres[i][[1]])
			trees.list.phy[i]<-read.simmap(text=trees.list[i])
		}
		
	####turn each tree into a phylo4d_ext class
	phy.ext.list<-vector("list", length(trees.list.phy))
		for(i in 1:length(trees.list.phy)){
			phy.ext.list[i]<-phyext(trees.list.phy[[i]])		}	
		
	####turn each tree into a brownie class	
	phy.brownie.list<-vector("list", length(phy.ext.list))
		for(i in 1:length(phy.ext.list)){
			phy.brownie.list[i]<-brownie(phy.ext.list[i])		}
		
	#####add the data to each tree
	phy.brownie.list.w.data<-vector("list", length(phy.brownie.list))
		for(i in 1:length(phy.brownie.list)){
			phy.brownie.list.w.data[i]<-addData(phy.brownie.list[i], tip.data=data, dataTypes=contData())
		}
		
		
	####now generate the "all" taxaset
	all_taxa<-grep("[a-z]", tipLabels(phy.brownie.list.w.data[[1]]), value=TRUE)  ##NOTE this is only for the junk tree -- need to fix this for the real species names 
		for(i in 1:length(phy.brownie.list.w.data)){
			taxasets(phy.brownie.list.w.data[[i]], taxnames="all")<-all_taxa		}

all.test.results1<-data.frame() #create an empty data frame as a repository of the final results
all.test.results2<-data.frame() #create an empty data frame as a repository of the final results

	
			
all.test.results1<-runNonCensored(phy.brownie.list.w.data, models=brownie.models.continuous()[2], treeloop=T, charloop=T)
all.test.results2<-runNonCensored(phy.brownie.list.w.data[[1]], models=brownie.models.continuous()[1], treeloop=T, taxset="all")



newRow<-data.frame(Tree=all.test.results2$Tree, Tree.weight =all.test.results2$"Tree weight", Tree.name=all.test.results2$"Tree name", Char=all.test.results2$Char, Model=all.test.results2$Model, LnL=all.test.results2$"-LnL", AIC=all.test.results2$AIC, AICc=all.test.results2$AICc, AncState=all.test.results2$AncState, Rate_in_state_0=all.test.results2$BMrate, Rate_in_state_1=all.test.results2$BMrate)

colnames(newRow)=c("Tree", "Tree weight", "Tree name", "Char", "Model", "-LnL", "AIC", "AICc", "AncState", "Rate_in_state_0", "Rate_in_state_1")

all.test.results<-rbind(all.test.results1, newRow)

		
dAIC=all.test.results$AIC-min(all.test.results$AIC)
dAICc=all.test.results$AICc-min(all.test.results$AICc)

AICweightraw=exp(-0.5*dAIC)
AICweight=AICweightraw/sum(AICweightraw)


AICcweightraw=exp(-0.5*dAICc)
AICcweight=AICcweightraw/sum(AICweightraw)

all.test.results<-cbind(all.test.results, dAIC, AICweight)
all.test.results<-cbind(all.test.results,dAICc, AICcweight)

colnames(all.test.results)[6]<-"-LnL"

return(all.test.results)

}


