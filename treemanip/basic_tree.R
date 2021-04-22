library(ape)
library(staphopia)
library(dplyr)
library(Biostrings)
library(seqinr)
library(phangorn)

TOKEN = "faccbc9aeb02dabcb6a1cc096e9e9a68c51d07c8"
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
core_index <- read.delim('nrd-gene-set.txt', sep = '\t', header=T)
core_index$id <- seq(1, dim(core_index)[1])
samples1 <- get_public_samples()
samples <- samples1$sample_id

#get 300 random samples from the staphopia database - first create random #s
#selectSamples <- sample(1:length(samples), 200, replace=F) #for random trees

#if using file read in (this gives sample ids, not indices, so it's sampleNums)
selectSamples <- read.delim('nrd-sample-set.txt')$sample_id #for nrd trees

#if using indices:
#sampleGenes <- get_variant_gene_sequence(as.numeric(samples[selectSamples]), annotation_ids = core_index$annotation_id)
#if we know which sample ids we want:
sampleGenes <- get_variant_gene_sequence(selectSamples,annotation_ids = core_index$annotation_id)

sampleGmr <- subset(sampleGenes, sample_id != 'reference')
gSamplemr <- sampleGmr %>% group_by(sample_id) %>% mutate(fullseq = paste0(sequence, collapse = ''))
gSamplemr <- gSamplemr[!duplicated(gSamplemr$sample_id),]
nexusDat <- list()
for(j in 1:length(selectSamples)){
    nexusDat[[j]] <- strsplit(tolower(gSamplemr$fullseq[j]), '')[[1]]
}
names(nexusDat) <- c(samples[selectSamples])
#nexus file is also a possibility for what to use - creating phydat/nexus/etc file type typically
#doesn't include names 

dna_list <- as.DNAbin.list(nexusDat)
phy_nexDat <- as.phyDat(dna_list)

#maybe we can just work in the dnabin format? 
#dna_dist <- dist.dna(dna_list, model = "JC69")
dna_dist <- dist.ml(phy_nexDat, model="JC69")

nexDat_UPGMA <- upgma(dna_dist)
nexDat_NJ  <- NJ(dna_dist)

plot(nexDat_UPGMA, main="UPGMA")

plot(nexDat_NJ, main = "Neighbor Joining")

#optimal stuff about optimal parsimony and pratchet
staph_optim <- optim.parsimony(nexDat_NJ, phy_nexDat)
staph_pratchet <- pratchet(phy_nexDat)

plot(staph_optim)
plot(staph_pratchet)
