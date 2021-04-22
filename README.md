# springrotation
R code from my spring rotation project with Dr. Daniel Weissman at Emory university to create phylogenetic trees using the Staphopia database.

The basic trees file can be used to pull from Staphopia various samples (either you can use sample IDs you know about or generated sample IDs - there are methods
for both, and they are marked in comments. Make sure to comment the one you're not using). The base block trees file is to be used after the basic trees file in 
order to generate the 3 Kb block trees mentioned in Nimwegen et. al. 2021. The other things are meant to be assistive for conversions and other things; the 
write.trees file includes some simple code to write trees to the current directory, while the published_samps file repeats the basic trees and base block trees
functions for published samples you might want to pull (just change the pubmed id). 
