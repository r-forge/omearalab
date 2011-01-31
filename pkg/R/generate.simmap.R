generate.simmap<-function (phy, taxa.vector, changeposition=0.5, digits = 10, tree.prefix = "") 
{
	#changeposition is where on the branch leading to the taxset to change from state 0 to state 1: at the beginning (changeposition=0, at the end (changeposition=1, at the midpoint (changeposition=0.5)))
    brl <- !is.null(phy$edge.length) #if TRUE, there are branch lengths
    nodelab <- !is.null(phy$node.label) #if TRUE, there are node labels
    phy$tip.label <- checkLabel(phy$tip.label) #gets rid of weird characters in species names
    if (nodelab) 
        phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "g", sep = "") #format to print out branch lengths (how many digits, always include the decimal point)
    cp <- function(x) {
        STRING[k] <<- x
        k <<- k + 1
    }
    getMRCA<-function(taxa.vector) {
    	mrcaNode<-which(phy$tip.label==taxa.vector[1])
    	for (taxonIndex in 1:length(taxa.vector)) {
    		mrcaNode=mrca(phy,full=TRUE)[mrcaNode,which(phy$tip.label==taxa.vector[taxonIndex])]
    	}	
    	return(mrcaNode)
    }
    #the getMRCA function will, for a given vector of taxon labels (and it currently assumes there are labels), will return the ape node number of the mrca of all those taxa
    
    #the next step is to make a vector of length equal to the number of internal and terminal nodes in the tree.
    #each entry of the vector will be one of three integers:
    #	0 = an edge (terminal or internal) not in the taxset
    #	1 = an edge completely in the taxset
    #	2 = an edge partly in the taxset, partly not in (so the stem branch leading to the MRCA of all the taxa in the taxset)
    
    #then, at the steps in the function that say "HERE IS ONE PLACE..." we want to instead do simmap.out(nodeNumber). simmap.out will format the simmap string and cp() it.
    #note that node number you pass to simmap.out should be somewhere from 1 to Ntip+Nnode. This differs from what seems to be done with add.internal, where it has to use the ind object. 
    
    #remember that simmap output: if the indicator vector is 0, it should be
    # {0,brlen}
    #if 1, it is {1,brlen}
    #if 2, it is {1,brlen*(1-changeposition) : 0,brlen*changeposition}
    #note that in the 2 case, you might want to do exceptions for changeposition==0 or ==1: cut out the part that has zero branch length
    simmap.out<-function(i) {
    	
    	
    	simmap.string<-paste("{",state.string,"}",sep="")
    	cp(simmap.string)	
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) 
                add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) 
                cp(",")
        }
        cp(")")
        if (nodelab && i > n) 
            cp(phy$node.label[i - n])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]])) #HERE IS ONE PLACE TO CHANGE FROM :2.646466 TO THE SIMMAP FORMAT. BUT WHICH STATE TO USE?
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i])) #HERE IS ONE PLACE TO CHANGE FROM :2.646466 TO THE SIMMAP FORMAT. BUT WHICH STATE TO USE?
        }
    }
    n <- length(phy$tip.label)
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent)) kids[[parent[i]]] <- c(kids[[parent[i]]], 
        children[i])
    ind <- match(1:max(phy$edge), phy$edge[, 2]) #I AM IND. BEWARE
    LS <- 4 * n + 5
    if (brl) 
        LS <- LS + 4 * n
    if (nodelab) 
        LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp(tree.prefix)
    cp("(")
    getRoot <- function(phy) phy$edge[, 1][!match(phy$edge[, 
        1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy)
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) 
            add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) 
            cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) 
            cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) 
            cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}
