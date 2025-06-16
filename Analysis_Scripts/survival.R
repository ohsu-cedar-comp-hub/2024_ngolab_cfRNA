library(survival)

# Use computed pdac scores and survival data frame to compute Kaplan-Meyer curves and p-values
survival_analysis = function(scores,survival_data,name,quan){
  survival = adjust_survival_score(survival_data,scores)
  
  #Create dataframes for survival plots 
  surv_data_c = make_survival_frame(survival,"CEDAR",quan,name)
  surv_data_b = make_survival_frame(survival,"BCC",quan,name)

  return(list(surv_data_c,surv_data_b))
}

adjust_survival_score = function(survival,scores){
  survival = survival[survival$SeqID %in% names(scores),]
  scores = scores[match(survival$SeqID,names(scores))]
  survival = cbind(survival,Discriminant.Score = scores)
  return(survival)
}


make_survival_frame = function(input_dataset,filter_source,quan,group){
  filtered_survival <- input_dataset[input_dataset$Cohort==filter_source,]
  surv_object <- Surv(time =filtered_survival$Followup.Duration..mon., event = filtered_survival$Outcome..NED..AWD..DOD== "DOD")
  Y <- surv_object
  
  quan2 = quantile(filtered_survival$Discriminant.Score,quan)
  kmfit <- survfit(Y~ (filtered_survival$Discriminant.Score>quan2))
  log_rank_test <- survdiff(formula = Y~ (filtered_survival$Discriminant.Score>quan2))
  
  calculated_p <- round(pchisq(log_rank_test$chisq, 1, lower=F), digits=3) 
  pval = c(paste0("p=",calculated_p),group)
  
  n1 = kmfit$strata[1]
  n2 = sum(kmfit$strata)
  surv_data = cbind(c(1,kmfit$surv[1:n1]),c(0,kmfit$time[1:n1]),paste0("Score < ",100*quan,"th percentile"),group )
  surv_data = rbind(surv_data,cbind(c(1,kmfit$surv[(n1+1):n2]),c(0,kmfit$time[(n1+1):n2]),paste0("Score > ",100*quan,"th percentile"),group))
  

  return(list(surv_data, pval, Y, filtered_survival))
}

  
