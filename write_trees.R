#paste0 removes unnecessary spaces
for (i in 1:201) {
  write.tree(njList[[i]],paste0("subset",i,".phy_phyml_tree.txt"))
}

#for UPGMA
for (i in 1:201) {
  write.tree(upgma_List[[i]],paste0("subset",i,".phy_phyml_tree.txt"))
}


#for pml
for (i in 1:length(optim_Trees)) {
  write.tree(optim_Trees[[i]],paste0("block",i,".phy_phyml_tree.txt"))
}