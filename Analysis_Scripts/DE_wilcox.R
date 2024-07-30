DE_wilcox = function(counts,groups,label_filt = NULL){
  if(!is.null(label_filt)){
    ids = groups %in% label_filt
    counts = counts[,ids]
    groups = groups[ids]
  }
  if(is.null(row.names(counts))){ row.names(counts) = 1:nrow(counts)}
  id1 = groups == groups[1]
  id2 = groups != groups[1]
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
