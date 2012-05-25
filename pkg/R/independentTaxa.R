#script to find independent species pairs: will return a data frame with pairs of taxa and the branch lengths for each
#none of these pairs share edges in the tree with any other pair
#Brian O'Meara, May 25, 2012

find.cherries<-function(phy) {
  focalNodes<-which(tabulate(phy$edge[which(phy$edge[,2]<=Ntip(phy)),1])==2)
  returnTaxa<-function(focalNode,phy) {
    return( c(phy$edge[which(phy$edge[,1]==focalNode ),2]))
  }
  returnBrlen<-function(focalTaxon,phy) {
    return( phy$edge.length[which(phy$edge[,2] == focalTaxon )] )
  }
  focalTaxa<-sapply(focalNodes,returnTaxa,phy)
  focalBrlen<-matrix(sapply(focalTaxa,returnBrlen,phy),ncol=2,byrow=TRUE)
  focalTaxa<-matrix(as.character(phy$tip.label[c(focalTaxa)]),ncol=2,byrow=TRUE)
  result<-data.frame(focalTaxa,focalBrlen,stringsAsFactors=FALSE)
  names(result)<-c("taxon1","taxon2","brlen1","brlen2")
  return(result)
}

independent.pairs<-function(phy,multi2di=TRUE) {
  pairs<-data.frame()
  phy<-multi2di(phy)
  while(Ntip(phy)>1) { #if odd number of taxa, will have one left over when all done
    cherries<-find.cherries(phy)
    pairs<-rbind(pairs,cherries)
    phy<-drop.tip(phy,tip=unlist(cherries[,1:2]))
  }
   return(pairs)
}