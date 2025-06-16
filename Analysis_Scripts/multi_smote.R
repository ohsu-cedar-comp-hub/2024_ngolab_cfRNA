#multiclass smote
#Interpolates data points using SMOTE to balance the 
#distribution of classes in input data.
#Inputs: 
#   X: Feature matrix (features x samples)
#   Y: Label vector
#Outputs: List of
#   X: Feature matrix with interpolated data points
#   Y: Label vector with interpolated labels
library(smotefamily)
multi_smote = function(X,Y){
  labels = unique(Y)
  ll = length(labels)
  max_class = names(which.max(table(Y)))
  #For each class, except for the majority class...
  for(c in 1:ll){
    if(labels[c]==max_class){
      next
    }
    #Perform 2 class SMOTE interpolation
    cids = Y == max_class | Y == labels[c]
    Xc = X[,cids]
    Yc = Y[cids]
    SM = SMOTE(data.frame(t(Xc)),Yc,K=min(5,sum(Y==labels[c])-1))
    
    #Add interpolated points to data
    Xc = SM$syn_data
    Yc = Xc[,ncol(Xc)]
    Xc = t(Xc[,-ncol(Xc)])
    row.names(Xc) = row.names(X)
    X = cbind(X,Xc)
    Y = c(Y,Yc)
  }
  return(list(X,Y))
}