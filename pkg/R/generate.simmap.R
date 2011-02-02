generate.simmap<-function (phy, taxa.vector, changeposition=0.5, digits = 10, tree.prefix = "") 
{
	#changeposition is where on the branch leading to the taxset to change from state 0 to state 1: at the beginning (changeposition=0, at the end (changeposition=1), at the midpoint (changeposition=0.5)))
    brl <- !is.null(phy$edge.length) #if TRUE, there are branch lengths
    nodelab <- !is.null(phy$node.label) #if TRUE, there are node labels
    phy$tip.label <- checkLabel(phy$tip.label) #gets rid of weird characters in tip names
    if (nodelab) 
        phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "g", sep = "") #format to print out branch lengths (how many digits, always include the decimal point)
    cp <- function(x) {
        STRING[k] <<- x
        k <<- k + 1
    }
    
            #the getMRCA function above will, for a given vector of taxon labels (and it currently assumes there are labels), return the ape node number of the mrca of all those taxa
    
    getMRCA<-function(phy,taxa.vector) {
    	mrcaNode<-which(phy$tip.label==taxa.vector[1])  #pic the label in the phylogeny that is the same as the first item in the taxa.vector and place it's index position into mrcaNode
    	for (taxonIndex in 1:length(taxa.vector)) {    #for each taxonIndex intance in taxa.vector
    		mrcaNode=mrca(phy,full=TRUE)[mrcaNode,which(phy$tip.label==taxa.vector[taxonIndex])]  #from the mrca(phy,full=T) result pull out the row that is equal to the mrcaNode position and any other tip taxon.label [mrcaNode, which(phy$ip.label==taxa.vector[taxonIndex])], for the column place of the mrca function
    	}	
    	return(mrcaNode)
    }
    
    
    
    	#the getImmediateDescendants, looks around the mrcaNode and pics out it's children
        getImmediateDescendants<-function(phy,mrcaNode){
    	   	vector.descendants<-phy$edge[which(phy$edge[,1]==mrcaNode),2]
       		return(vector.descendants)		
       }
       
              
       #the getAllDescendants uses the results of the getImmediateDescendants() 
       getAllDescendants<-function(phy,mrcaNode){
       		vector.alldescendants<-getImmediateDescendants(phy,mrcaNode) #first get results 
       		to.examine=vector.alldescendants   #stick the getImm.Dec() res into a new vector
       		while(length(to.examine)>0) {  #while the immediate descendants are greater than 0
       			vector.nextdescendants=getImmediateDescendants(phy,to.examine[1]) #get Immediate desc for the first element in to.examine stick it in vector.nextdesc.
       			if (length(vector.nextdescendants)>0 ) {   #if vector.nextdesc. continues to have more than zero items in it
       				to.examine<-append(to.examine,vector.nextdescendants) #then append to.examine, with the results of vector.nextdesc, stick it into to.examine
       				vector.alldescendants<-append(vector.alldescendants,vector.nextdescendants) #append vector.alldesc, with the results of vector.nextdesc, into vector all.desc
       			}	
       			if (length(to.examine)>1) {   #if to.examine is greater than one
	       			to.examine<-to.examine[2:length(to.examine)]   #take the second element and on of to.examine and stick it in to.examine, i.e. drop the first element of to.examine
       			}
       			else {
       				to.examine=c()	#if to.examine is not greater than one, return an empty to.examine vector, i.e clear it out. 
       			}
       		}
       		return(vector.alldescendants)
       	
       	}
        
    
        #the next step is to make a vector of length equal to the number of internal and terminal nodes in the tree.
    #each entry of the vector will be one of three integers:
    #	0 = an edge (terminal or internal) not in the taxset
    #	1 = an edge completely in the taxset
    #	2 = an edge partly in the taxset, partly not in (so the stem branch leading to the MRCA of all the taxa in the taxset)
    
    nodes.assignment<-rep(0, Nnode(phy) + Ntip(phy))
    
    nodes.assignment[mrcaNode]=2      #assign mrcaNode element to be 2 to show where the break is
    
    nodes.assignment[getAllDescendants(phy, mrcaNode)]=1   #this might be out of the scope for the getAllDescendants function above where vector.alldescendatns gets created   
    

    
    
    #then, at the steps in the function that say "HERE IS ONE PLACE..." we want to instead do simmap.out(nodeNumber). simmap.out will format the simmap string and cp() it.
    #note that nodeNumber you pass to simmap.out should be somewhere from 1 to Ntip+Nnode. This differs from what seems to be done with add.internal, where it has to use the ind object. 
    
    #remember that simmap output: 
    # if the indicator vector is 0, it should be {0,brlen}
    #if 1, it is {1,brlen}
    #if 2, it is {1,brlen*(1-changeposition) : 0,brlen*changeposition}
    #note that in the 2 case, you might want to do exceptions for changeposition==0 or ==1: cut out the part that has zero branch length
    
    
    
    simmap.out<-function(nodes.assignment, changeposition, i) {   # simmap.out will format the simmap string and cp() it. need to feed it a nodeNumber
    	state.type<-nodes.assignment[i]
    	if (state.type==0) {
    		cp(paste("{0,", sprintf(f.d, phy$edge.length[i]), "}", sep=""))
    		}
    		else if (state.type==1){ 
    			cp(paste("{1," sprintf(f.d, phy$edge.length[i]), "}", sep=""))
    			}
    		else if	(state.type==2){ 
    			cp(paste("{1," sprintf(f.d, phy$edge.length[i]*(1-changeposition)),":","0,", sprintf(f.d, phy$edge.length[i]*(changeposition)), "}", sep=""))
    			
    			
    			#modify  this last cp statement for if changeposition =0 or =1
    			
    			
    			}
    			
    	#remember the colon is already output by cp(":")
    	 #going to be callig one node at a time
    	  
    	  
    		
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
            simmap.out(i) #HERE IS ONE PLACE TO CHANGE FROM :2.646466 TO THE SIMMAP FORMAT. BUT WHICH STATE TO USE?
        }
    }
    
    
    
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i])) #HERE IS ONE PLACE TO CHANGE FROM :2.646466 TO THE SIMMAP FORMAT. BUT WHICH STATE TO USE?
        }
    }
    
    
    
    ####tree traversal to assign children of each node, parent nodes, and tip taxa identification
    
    n <- length(phy$tip.label)
    
    parent <- phy$edge[, 1]
    
    children <- phy$edge[, 2]
    
    kids <- vector("list", n + phy$Nnode)
    
    for (i in 1:length(parent)) kids[[parent[i]]] <- c(kids[[parent[i]]], 
        children[i])
       
    ####tree traversal to assign children of each node, parent nodes, and tip taxa identification
     
        
        
        
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
