### Decompose matrix 'Y' into 'f' factors.
### Runs 'reps' random restarts, returns best performing iteration.
###     Performance scored by reconstruction loss of input matrix.
intrinsic_factors = function(Y,f,reps){
  library(MASS)
  n = dim(Y)[1]
  p = dim(Y)[2]
  max_iter = 50
  best_recon = 1000
  Wb = 0
  Hb = 0
  for(rep in 1:reps){
    W = matrix(runif(n*f),ncol=f)
    W = 1000000*t(t(W)/colSums(W))
    H = matrix(0.5,ncol=p,nrow=f)
    conv = FALSE
    Wn = W
    Hn = H
    iter = 0
    while(!conv){
      iter = iter+1
      Wn = W*(Y%*%t(H))/(W%*%H%*%t(H))
      Wn[is.nan(Wn)] = 0
      for(j in 1:p){
        Hn[-1,j] = H[-1,j]*(t(Wn)%*%Y)[-1,j]/((t(Wn)%*%Wn%*%H)[-1,j]) 
        if(sum(Hn[-1,j]) > .999){
          Hn[-1,j] = .999*Hn[-1,j]/sum(Hn[-1,j])
        }
        
        Hn[1,j] = 1-sum(Hn[-1,j])
      }
      
      Hn[is.nan(Hn)] = 0
      diff_H = sum(abs(H-Hn))/sum(H)
      diff_W = sum(abs(W-Wn))/sum(W)
      if(diff_W+diff_H < .002 || iter > max_iter){
        conv=TRUE
      }
      W = Wn
      H = Hn
    }
    recon = norm(Y-W %*% H,type="F")/norm(Y,type="F")
    if(recon < best_recon){
      best_recon = recon
      Wb = W
      Hb = H
    }
  }
  R = list(Wb,Hb,best_recon)
  return(R)
}

### Adjust factors 'fixed' to data 'Y'
### Adds new entries for each feature in 'Y' not in 'fixed'
### Only computes factor values for new features.
### Runs 'reps' random restarts, returns best performing iteration.
###     Performance scored by reconstruction loss of input matrix.
adjust_factors = function(Y,fixed,reps){
  library(MASS)
  f = ncol(fixed)
  n = dim(Y)[1]
  p = dim(Y)[2]
  fixed = fixed[row.names(fixed) %in% row.names(Y),]
  fids = match(row.names(fixed),row.names(Y))
  fixed = fixed * rowMeans(Y)[fids]
  max_iter = 50
  best_recon = 1000
  Wb = 0
  Hb = 0
  for(rep in 1:reps){
    W = matrix(runif(n*f),ncol=f)
    W = 1000000*t(t(W)/colSums(W))
    W[fids,] = fixed
    H = matrix(0.5,ncol=p,nrow=f)
    conv = FALSE
    Wn = W
    Hn = H
    iter = 0
    while(!conv){
      iter = iter+1
      Wn = W*(Y%*%t(H))/(W%*%H%*%t(H))
      Wn[fids,] = fixed
      Wn[is.nan(Wn)] = 0
      for(j in 1:p){
        Hn[-1,j] = H[-1,j]*(t(Wn)%*%Y)[-1,j]/((t(Wn)%*%Wn%*%H)[-1,j]) 
        if(sum(Hn[-1,j]) > .999){
          Hn[-1,j] = .999*Hn[-1,j]/sum(Hn[-1,j])
        }
        
        Hn[1,j] = 1-sum(Hn[-1,j])
      }
      
      Hn[is.nan(Hn)] = 0
      diff_H = sum(abs(H-Hn))/sum(H)
      diff_W = sum(abs(W-Wn))/sum(W)
      if(diff_W+diff_H < .002 || iter > max_iter){
        conv=TRUE
      }
      W = Wn
      H = Hn
    }
    recon = norm(Y-W %*% H,type="F")/norm(Y,type="F")
    if(recon < best_recon){
      best_recon = recon
      Wb = W
      Hb = H
    }
  }
  R = list(Wb,Hb,best_recon)
  return(R)
}

