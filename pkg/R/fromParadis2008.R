#finding a clade (Paradis, 2008: FormatTreeR)
#in previous example Paradis was walking up a tree from the root node

For an arbitrary node, say nod, the approach above may be adapted to the
following algorithm:
1. Find the node ancestor of nod, store its number in anc, and store the
number of the branch in i.
2. Find the next occurrence of anc in tr$edge[, 1], store it in j.

The clade descending from nod is given by the rows i+1 to j−1 of tr$edge.
This algorithm coded in R is:

i <- which(phy$edge[, 2] == nod)  # find position (row) of nod in column 2
anc <- phy$edge[i, 1]    #what is the value in column 1, given row i; gives you ancestor of nod
tmp <- which(phy$edge[, 1] == anc) #which rows (i.e ancestors, col 1) in phy are equal to anc (i.e what are the total nodes that are subtended by anc)
j <- tmp[which(tmp == i) + 1] #in tmp which item matches i? then add 1 to it. then then store it in j; gives the row i +1 that is an anc for nod
phy$edge[(i+1):(j-1), ] #now need to output everything in between the phy rows i + 1, and j - 1, which are the branches for the clade which descends from nod.