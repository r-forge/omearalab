library(ape)
library(geiger)
#transforms trees so substitution rate evolves under brownian motion. 
#The state at each node is used to transform the branch length of the input tree.
heritableTransform<-function(phy,rate=1) {
  tips<-abs(1+sim.char(phy,matrix(rate),model="brownian",nsims=1))
  ancstates<-(ace(tips,phy))$ace
  ancstates<-ancstates[-1] #remove root value
  matches<-match(as.numeric(names(ancstates)),phy$edge[,2])
  phy$edge.length[matches]<-phy$edge.length[matches]*ancstates
  tip.id<-match(row.names(tips),phy$tip.label)
  matches<-match(tip.id,phy$edge[,2])
  phy$edge.length[matches]<-phy$edge.length[matches]*tips
  return(phy)
}
