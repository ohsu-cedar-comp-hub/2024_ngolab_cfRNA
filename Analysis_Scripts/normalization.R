# Counts per million normalization of counts matrix (by column)
cpm = function(C){
  C = t(t(C)/colSums(C))*1e6
  return(C)
}

# Trimmed mean of M-value normalization of counts matrix (by column)
tmm = function(C){
  library(edgeR)
  scale.factors = calcNormFactors(C,lib.size=NULL, method = "TMM")
  C = t(t(C)/(scale.factors))
  return(C)
}

# Our cf-normalization method
# Takes in a counts matrix and set of gene factors
#     Factors should have an Intrinsic and Extrinsic component
cf_norm = function(C,intrinsic_genes){
  intrinsic_genes = intrinsic_genes[intrinsic_genes[,1] %in% row.names(C),]
  C = C[intrinsic_genes[,1],]
  Ci = C * intrinsic_genes$Intrinsic
  Cb = C * intrinsic_genes$Extrinsic
  
  fact_i = colSums(Ci)
  fact_i = fact_i/median(fact_i)
  fact_b = colSums(Cb)
  fact_b = fact_b/median(fact_b)
  
  FM = as.matrix(intrinsic_genes[,c(2,3)]) %*% as.matrix(rbind(1/fact_i,1/fact_b))
  
  C = C * FM
  return(C)
}

# Transcripts per million normalization of counts matrix
# Assumes the file g_to_e.csv, which has the transcript length values, is in the working directory
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
