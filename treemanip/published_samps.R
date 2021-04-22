#code for using published samples from staphopia
TOKEN = "faccbc9aeb02dabcb6a1cc096e9e9a68c51d07c8"
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
core_index <- read.delim('nrd-gene-set.txt', sep = '\t', header=T)
core_index$id <- seq(1, dim(core_index)[1])
published <- get_published_samples()
publishedids <- published$sample_id
publications <- get_publications(publishedids)

#got publications somehow with pmids - forgot how, so check with rohan or 
#remember

#indices in publications don't match indices in published. So i need to go off of sample ids ->
licationsIDs <- publications$sample_id[which(publications$pmid == 22393007)]#grab sample ids with desired pmid from publications list

lishedIDs <- published$sample_id #grab sample ids from the published list
indexList <- which(is.element(lishedIDs,licationsIDs)) #get indices of samples with matching sample ids from published

#now collect the gene sequences we want and do the standard organizational steps.
groupGenes <- get_variant_gene_sequence(as.numeric(lishedIDs[indexList]), annotation_ids = core_index$annotation_id)

groupGmr <- subset(groupGenes, sample_id != 'reference')
gGroupmr <- groupGmr %>% group_by(sample_id) %>% mutate(fullseq = paste0(sequence, collapse = ''))
gGroupmr <- gGroupmr[!duplicated(gGroupmr$sample_id),]
groupNexusDat <- list()
for(j in 1:length(indexList)){
  groupNexusDat[[j]] <- strsplit(tolower(gGroupmr$fullseq[j]), '')[[1]]
}
names(groupNexusDat) <- c(published[indexList]$sample_id)

#now using phydat to convert to phangorn
phy_gNexDat <- phyDat(groupNexusDat, type = "DNA", levels = NULL)

gDna_dist <- dist.ml(phy_gNexDat, model="JC69")

gnexDat_UPGMA <- upgma(gDna_dist)
gnexDat_NJ  <- NJ(gDna_dist)

plot(gnexDat_UPGMA, main="UPGMA")

plot(gnexDat_NJ, main = "Neighbor Joining")


#now doing base blocks!


#starting from the point where nexusdat has already been created for 100 samples (see basic_tree)
gbigList <- list()
gno <- -2999
gnew <- 0
#create list of of lists of samples for different segments of the core genome
for (j in 1:200){
  gno <- gno + 3000
  gnew <- gnew + 3000
  gnewList <- list()
  for (i in 1:length(groupNexusDat)){
    gnewList[[i]] <- groupNexusDat[[i]][gno:gnew]
  }
  gbigList[[j]] <- gnewList
  names(gbigList[[j]]) <- names(groupNexusDat)
}
#convert to phydat for use with phangorn
for(i in 1:length(gbigList)) {
  gbigList[[i]] <- phyDat(gbigList[[i]], type = "DNA", levels = NULL)
}


#now create a list of dna distance matrices
gdistList <- list()
for (i in 1:length(gbigList)) {
  gdistList[[i]] <- dist.ml(gbigList[[i]], model = "JC69")
}


#then combine with all the trees - should end up as type "multiphylo". Dropping
#negative edge lengths using drop labels. 
gnj_List <- c(gnexDat_NJ)
for (i in 1:length(gdistList)){
  gnj_List <- c(gnj_List, NJ(gdistList[[i]]))
}

gupgma_List <- c(gnexDat_UPGMA)
for (i in 1:length(gdistList)) {
  gupgma_List <- c(gupgma_List, upgma(gdistList[[i]]))
}

#creating optimal trees?
goptim_Trees <- c(pml(gnexDat_NJ, phy_gNexDat)$tree)
for (i in 2:length(gnj_List)) {
  goptim_Trees <- c(goptim_Trees, pml(gnj_List[[i]],gbigList[[i-1]])$tree)
}

#finally, create a distance matrix, using the total core genome as first entry
dist.topo(gupgma_List, method = "score") #for viewing
dmatElife <- dist.topo(gupgma_List, method = "score") #for later analysis

#when testing subtrees
gSubs <- subtrees(goptim_Trees[[1]])
gSubs[[80]]$tip.label
selectSamples <- as.numeric(gSubs[[80]]$tip.label)
