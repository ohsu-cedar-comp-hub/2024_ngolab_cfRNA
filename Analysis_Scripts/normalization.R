cpm = function(C){
  C = t(t(C)/colSums(C))*1e6
  return(C)
}

tmm = function(C){
  library(edgeR)
  scale.factors = calcNormFactors(C,lib.size=NULL, method = "TMM")
  C = t(t(C)/(scale.factors))
  return(C)
}

tpm = function(D){
  map = read.csv("g_to_e.csv")
  genes = row.names(D)
  #Get matching rows
  D = D[genes %in% map[,2],]
  ids = match(genes,map[,2])
  map = map[ids,]
  #divide counts by gene lengths
  D = D/map[,5]
  D = D[map[,5]>0,]
  #throw into CPM
  C = cpm(D)
  
  return(C)
}

intrinsic_norm = function(C,intrinsic_genes){
  ids = row.names(C) %in% intrinsic_genes
  fact = colSums(C[ids,])
  C = t(t(C)/fact)
  C = cpm(C)
  return(C)
}