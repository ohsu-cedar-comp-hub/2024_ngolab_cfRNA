library(e1071)
svm_decon = function(C,S){
  res = matrix(0,nrow=dim(C)[2],ncol=dim(S)[2])
  row.names(res) = colnames(C)
  colnames(res) = colnames(S)

  C = scale(C)
  C[is.nan(C)]=0
  S = scale(S)
  
  for(i in 1:dim(C)[2]){
    b_val = Inf
    for(v in c(0.05, 0.1, 0.15, 0.25, 0.5, 0.75)){
      model = svm(C[,i] ~ S,type="nu-regression",kernal="linear",nu=v)
      #get coefficients
      coef <- t(model$coefs) %*% model$SV
      #remove negative coefficients, sum to 1
      coef[coef<0]=0
      coef = coef/sum(coef)
      #SSE
      resid = sum((C[,i]-S%*%t(coef))^2)
      if(resid < b_val){
        b_val = resid
        b_coef = coef
      }
    }
    res[i,] = b_coef
  }
  return(res)
}

