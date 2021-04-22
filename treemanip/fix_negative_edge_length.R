fix_negative_edge_length <- function(nj.tree) {
  edge_infos <- as.data.table(cbind(nj.tree$edge, nj.tree$edge.length)) 
  colnames(edge_infos) <- c('from', 'to', 'length')
  nega_froms <- edge_infos[length < 0, sort(unique(from))]
  for (nega_from in nega_froms) {
    minus_length <- edge_infos[from == nega_from, ][order(length)][1, length]
    edge_infos[from == nega_from, length := length - minus_length]
    edge_infos[to == nega_from, length := length + minus_length]
  }
  nj.tree$edge.length <- edge_infos$length
  nj.tree
}
#Nam, Le. (2020). Re: How to correct negative branches from neighbor-joining method?.
#Retrieved from: 
#https://www.researchgate.net/post/How-to-correct-negative-branches-from-neighbor-joining-method/5fa2790926cfd7752128b32d/citation/download. 
