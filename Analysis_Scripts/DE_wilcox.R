# Wilcoxon rank sum test for differential genes
# Takes count table (matrix) and binary labels as inputs
# By default labels are assumed to be boolean, with the values of TRUE treated as the first class
# A character vector of labels can also be passed
#     In this case comp_label should be the string of the positive label, with all others treated as negative
DE_wilcox = function(counts,groups,label_filt = NULL,comp_label=TRUE){
  if(!is.null(label_filt)){
    ids = groups %in% label_filt
    counts = counts[,ids]
    groups = groups[ids]
  }
  if(is.null(row.names(counts))){ row.names(counts) = 1:nrow(counts)}
  id1 = groups == comp_label
  id2 = groups != comp_label
  D = matrix(0,nrow=nrow(counts),ncol=4)
  #row.names(D) = row.names(counts)
  
  for(i in 1:nrow(D)){
    x = as.numeric(counts[i,])
    p = wilcox.test(x[id1],x[id2])$p.value
    if(is.na(p)){p=1}
    change = log2(mean(x[id1]+1))-log2(mean(x[id2]+1))
    if(is.na(change)){
      change = 0
    }
    D[i,] = c(row.names(counts)[i],change,mean(x),p)
  }
  D = cbind(D,p.adjust(D[,4],method="fdr"))
  colnames(D) = c("genes","LogFC","LogCPM","PValue","FDR")
  D = data.frame(D)
  return(D)
}
