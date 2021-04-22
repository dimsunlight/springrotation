#starting from the point where nexusdat has already been created for a given number of samples (see basic_tree)
bigList <- list()
no <- -2999
new <- 0
noList <- c(no)
newList <- c(new)
#create list of of lists of samples for different segments of the core genome
for (j in 1:200){
  no <- no + 3000
  new <- new + 3000
  noList <- c(noList,no)
  newList <- c(newList,no)
  n1List <- list()
  #convert to phydat within list creation to avoid bugs from converting later
  for (i in 1:length(nexusDat)){
    n1List[[i]] <- nexusDat[[i]][no:new]
  }
  bigList[[j]] <- n1List
  names(bigList[[j]]) <- names(nexusDat)
}

#convert to nexus dat
for (i in 1:length(bigList)){
  bigList[[i]] <- phyDat(bigList[[i]], type="DNA", levels = NULL)
}

#now create a list of dna distance matrices
distList <- list()
for (i in 1:length(bigList)) {
  distList[[i]] <- dist.ml(bigList[[i]], model = "JC69")
}


#then combine with all the trees - should end up as type "multiphylo. Dropping
#negative edge lengths using drop labels. 
nj_List <- c(nexDat_NJ)
for (i in 1:length(distList)){
  nj_List <- c(nj_List, NJ(distList[[i]]))
}

upgma_List <- c(nexDat_UPGMA)
for (i in 1:length(distList)) {
  upgma_List <- c(upgma_List, upgma(distList[[i]]))
}

#creating optimal trees?
optim_Trees <- c(pml(nexDat_NJ, phy_nexDat)$tree)
for (i in 2:length(nj_List)) {
  optim_Trees <- c(optim_Trees, pml(nj_List[[i]],bigList[[i-1]])$tree)
}

#finally, create a distance matrix, using the total core genome as first entry
dist.topo(upgma_List, method = "score") #for viewing
dmatElife <- dist.topo(upgma_List, method = "score") #for later analysis
library(profvis)
for (i in 1:length(optim_Trees)){
  plot(optim_Trees[[i]]) 
  pause(.5)
  }