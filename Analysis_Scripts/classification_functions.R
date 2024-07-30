classifier_suite = function(counts_train,labels_train,meta_train,counts_val,labels_val,meta_val,params=list()){
  BvsP_train = meta_train$Group %in% c("Benign","PDAC")     #Benign and PDAC samples
  BvsP_val = meta_val$Group %in% c("Benign","PDAC")
  NCvsP_train = meta_train$Group != "Other Cancer"          #PDAC and other non-cancer samples
  NCvsP_val = meta_val$Group != "Other Cancer"
  
  #10-fold cross validation
  kfold_b_vs_p = k_fold_validation(counts_train[,BvsP_train],labels_train[BvsP_train],10,alg="RF",params=params)
  kfold_nc_vs_p = k_fold_validation(counts_train[,NCvsP_train],labels_train[NCvsP_train],10,alg="RF",params=params)
  kfold_all = k_fold_validation(counts_train,labels_train,10,alg="RF",params=params)
  
  #External validation
  val_b_vs_p = test_validation(counts_train[,BvsP_train],
                               labels_train[BvsP_train],
                               counts_val[,BvsP_val],
                               labels_val[BvsP_val],
                               alg="RF",params=params)
  val_nc_vs_p = test_validation(counts_train[,NCvsP_train],
                                labels_train[NCvsP_train],
                                counts_val[,NCvsP_val],
                                labels_val[NCvsP_val],
                                alg="RF",params=params)
  val_all = test_validation(counts_train,
                            labels_train,
                            counts_val,
                            labels_val,
                            alg="RF",params=params)
  
  return(list(kfold_b_vs_p,kfold_nc_vs_p,kfold_all,val_b_vs_p,val_nc_vs_p,val_all))
}



#Perform cross-validation of data with different choices of parameters
#parameter lists are defined for LDA, RF, SVM
#Returns average per-class AUC for each configuration of parameters
param_search = function(X,Y,k,alg,ratio=0,useSmote=TRUE){
  results = c()
  if(alg=="LDA"){
    method = c("mle","moment","mve","t")
    nu = c(3,5,10)
    for(m in method){
      for(n in nu){
        auc = k_fold_validation(X,Y,k,ratio=ratio,alg=alg,useSmote=useSmote,params=list(method=m,nu=n))[[1]]
        results = rbind(results,c(mean(auc),m,n))
      }
    }
    colnames(results) = c("avg_auc","method","nu")
  } else if(alg=="SVM"){
    cost = c(.1,.5,1,5,10)
    for(c in cost){
      auc = k_fold_validation(X,Y,k,ratio=ratio,alg=alg,useSmote=useSmote,params=list(cost=c))[[1]]
      results = rbind(results,c(mean(auc),c))
    }
    colnames(results) = c("avg_auc","cost")
  } else if(alg=="RF"){
    mtry = c(1,5,10)
    ntrees = c(250,500,1000)
    nodesize = c(1,5,10)
    replace = c(TRUE,FALSE)
    for(m in mtry){
      for(nt in ntrees){
        for(ns in nodesize){
          for(r in replace){
            params = list(mtry=m,replace=r,ntrees=nt,nodesize=ns)
            auc = k_fold_validation(X,Y,k,ratio=ratio,alg=alg,useSmote=useSmote,params=params)[[1]]
            results = rbind(results,c(m,r,ns,nt,mean(auc)))
          }
        }
      }
    }
    colnames(results) = c("mtry","replace","nodesize","ntrees","avg_auc")
  }
  return(results)
}


k_fold_validation = function(X,Y,k,ratio=0,alg="RF",useSmote=TRUE,useLog=TRUE,params=list(),labels = unique(Y)){
  library(randomForest)
  library(e1071)
  library(MASS)
  library(pROC)
  source("multi_smote.R")
  
  ids = Y %in% labels
  Y = Y[ids]
  X = X[,ids]
  
  if(useLog){
    X = log(X+1)
  }
  
  set.seed(5)
  if(ratio==0){
    #k-fold splitting
    mult = ceiling(length(Y)/k)
    fold = rep(0,length(Y))
    for(label in labels){
      fold[Y==label] = sample(rep(1:k,mult)[1:sum(Y==label)])
    }
  }
  avg_auc = matrix(0,nrow=k,ncol=length(labels))
  colnames(avg_auc) = labels
  
  pred = matrix(0,nrow=length(Y),ncol=length(labels))
  colnames(pred) = labels
  row.names(pred) = colnames(X)
  n_pred = rep(0,length(Y))
  
  all_feats = c()
  #for each fold
  for(i in 1:k){
    #split training and validation, by groups
    if(ratio != 0){
      set.seed(i*5)
      fold = rep(0,length(Y))
      fold[sample(1:length(Y),round((1-ratio)*length(Y)),replace = FALSE)] = i
    }
    
    t_ids = c()
    v_ids = which(fold==i)
    t_ids = which(fold!=i)

    X_train = X[,t_ids]
    X_val = X[,v_ids]
    Y_train = Y[t_ids]
    Y_val = Y[v_ids]
    
    if(useSmote){
      sm_res = multi_smote(X_train,Y_train)
      X_train = sm_res[[1]]
      Y_train = sm_res[[2]]
    }
    
    #Fit model and predict
    if(alg == "LDA"){
      #mod = lda(t(X_train),factor(Y_train))
      mod = do.call(lda,c(list(x=t(X_train),grouping=factor(Y_train)),params))
      PM = predict(mod,t(X_val))$posterior
    } else if(alg == "SVM"){
      #mod = svm(t(X_train),factor(Y_train),kernal="linear",probability=TRUE)
      mod = do.call(svm,c(list(x=t(X_train),y=factor(Y_train),kernal="linear",probability=TRUE),params))
      PM = attr(predict(mod,t(X_val),probability=TRUE),"probabilities")
    } else if(alg == "RF"){
      mod = do.call(randomForest,c(list(x=t(X_train),y=factor(Y_train)),params))
      #mod = randomForest(t(X_train),factor(Y_train))
      PM = predict(mod,t(X_val),type="prob")
      
    }
    
    pred[v_ids,] = pred[v_ids,]+PM[,labels]
    n_pred[v_ids] = n_pred[v_ids]+1
    for(j in 1:length(labels)){
      if(length(unique(Y_val))>1){
        avg_auc[i,j] = auc(roc(factor(Y_val==labels[j]),PM[,labels[j]]))
      }
    }
  }
  
  for(i in 1:ncol(pred)){
    pred[,i] = pred[,i]/n_pred
  }
  
  #evaluate
  AUC = rep(0,length(labels))
  names(AUC) = labels
  for(i in 1:length(labels)){
    AUC[i] = auc(roc(factor(Y==labels[i]),pred[,i]))
  }
  
  if(ratio==0){
    return(list(AUC,pred,all_feats))
  } else{
    return(list(colMeans(avg_auc),pred,all_feats))
  }
}

test_validation = function(X_train,Y_train,X_val,Y_val,alg="RF",useSmote=TRUE,useLog=TRUE,params=list(),labels=unique(Y_train)){
  library(randomForest)
  library(e1071)
  library(MASS)
  library(pROC)
  source("multi_smote.R")
  set.seed(50)
  
  ids1 = Y_train %in% labels
  ids2 = Y_val %in% labels
  Y_train = Y_train[ids1]
  X_train = X_train[,ids1]
  Y_val = Y_val[ids2]
  X_val = X_val[,ids2]
  
  if(useLog){
    X_train = log(X_train+1)
    X_val = log(X_val+1)
  }
  
  pred = matrix(0,nrow=length(Y_val),ncol=length(labels))
  colnames(pred) = labels
  row.names(pred) = colnames(X_val)
  
  
  X_train_o = X_train
  
  #imputation
  if(useSmote){
    sm_res = multi_smote(X_train,Y_train)
    X_train = sm_res[[1]]
    Y_train = sm_res[[2]]
  }
  
  #Fit model and predict
  if(alg == "LDA"){
    mod = do.call(lda,c(list(x=t(X_train),grouping=factor(Y_train)),params))
    PM = predict(mod,t(X_val))$posterior
    PMt = predict(mod,t(X_train_o))$posterior
  } else if(alg == "SVM"){
    mod = do.call(svm,c(list(x=t(X_train),y=factor(Y_train),kernal="linear",probability=TRUE),params))
    PM = attr(predict(mod,t(X_val),probability=TRUE),"probabilities")
    PMt = attr(predict(mod,t(X_train_o),probability=TRUE),"probabilities")
  } else if(alg == "RF"){
    mod = do.call(randomForest,c(list(x=t(X_train),y=factor(Y_train)),params,importance=TRUE))
    PM = predict(mod,t(X_val),type="prob")
    PMt = predict(mod,t(X_train_o),type="prob")
  }
  
  pred = PM[,labels]
  predt = PMt[,labels]
  
  #evaluate
  AUC = rep(0,length(labels))
  names(AUC) = labels
  if(length(unique(Y_val))>1){
    for(i in 1:length(labels)){
      AUC[i] = auc(roc(factor(Y_val==labels[i]),pred[,i]))
    }
  }
  return(list(AUC,pred,predt))
}
