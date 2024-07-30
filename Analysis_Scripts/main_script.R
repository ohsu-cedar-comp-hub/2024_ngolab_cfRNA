#### Main script to reproduce PDAC classification analysis ####
source("normalization.R")
source("DE_wilcox.R")
source("deconvolve_cfRNA.R")
source("paper_figures.R")
source("classification_functions.R")
source("survival.R")

#### Load data and meta data ####
rna_counts = read.csv("pdac_genecount.csv",row.names=1)
meta = read.csv("pdac_meta.csv")

#Combine islet cell tumor diagnosis with cancer_other
meta$Group[meta$Group=="Islet Cell Tumor"] = "Other Cancer"
meta$Group[meta$Group=="Cancer_other"] = "Other Cancer"
meta$Group[meta$Group=="Benign Pancreas"] = "Benign"
meta$Source[meta$Source=="BCC_2018"] = "BCC_2019"

#Set color schemes
group_colors = c("forestgreen","blue","gold","purple","red")
group_levels = c("Benign","Pancreatitis","IPMN","PDAC","Other Cancer")
names(group_colors) = group_levels

#prep id sets for subset analysis
ids_B_vs_P = meta$Group %in% c("Benign","PDAC")                 #Benign and PDAC samples
ids_NC_vs_P = meta$Group != "Other Cancer"                      #PDAC and other non-cancer samples
ids_all = rep(TRUE,nrow(meta))                                  #All samples
ids_noLivMet = meta$Met_Site != "liver" | meta$Group!="PDAC"    #All samples without PDAC - liver metastasis
ids_stage12 = !meta$Stage %in% c(3,4) | meta$Group!="PDAC"      #All samples without stage 3-4 PDAC
ids_stage34 = !meta$Stage %in% c(1,2) | meta$Group!="PDAC"      #All samples without stage 1-2 PDAC
ids_male = meta$Gender.Flag == "M"
ids_ca19 = !is.na(meta$CA19_9)                                  #Samples that have CA 19-9 biomarker measurements
ids_train = meta$Source=="CEDAR_2020"
ids_val = meta$Source == "BCC_2019"

#### Perform intrinsic gene normalization + TMM ####
intrinsic_genes = read.csv("Intrinsic_score.csv")
norm_gene_set = intrinsic_genes[1:100,1]
rna_counts_norm = tmm(intrinsic_norm(rna_counts,norm_gene_set))

### Separate into training and validation sets ####
train_counts_norm = rna_counts_norm[,ids_train]
val_counts_norm = rna_counts_norm[,ids_val]
train_labels = meta$Group[ids_train]
val_labels = meta$Group[ids_val]

#remove low variance genes, computed on training set
V = apply(train_counts_norm,1,var)
train_counts_norm = train_counts_norm[V > 1,]
val_counts_norm = val_counts_norm[V > 1,]

#### Feature selection and filtering ####
DE = DE_wilcox(train_counts_norm,train_labels=="PDAC")
DE_genes = DE[DE[,"FDR"]<.05,1]


#### Deconvolution analysis of selected geneset ####
atlas = read.csv("Tissue_RNA_Atlas.csv",row.names=1)
#set minimum expression to 1, for deconvolution stability
atlas[atlas < 1] = 1
genes = DE_genes[DE_genes %in% row.names(atlas)]
atlas = atlas[genes,]
dc_svm_all = svm_decon(cbind(train_counts_norm[genes,],val_counts_norm[genes,]),atlas)
dc_svm_all = dc_svm_all / rowSums(dc_svm_all)
dc_svm_noLivMet = svm_decon(cbind(train_counts_norm[genes,ids_noLivMet[ids_train]],val_counts_norm[genes,ids_noLivMet[ids_val]]),atlas)
dc_svm_noLivMet = dc_svm_noLivMet / rowSums(dc_svm_noLivMet)
dc_svm_stage12 = svm_decon(cbind(train_counts_norm[genes,ids_stage12[ids_train]],val_counts_norm[genes,ids_stage12[ids_val]]),atlas)
dc_svm_stage12 = dc_svm_stage12 / rowSums(dc_svm_stage12)
dc_list = list(dc_svm_all,dc_svm_noLivMet,dc_svm_stage12)

#### Classification ####
train_labels_oneClass = train_labels
train_labels_oneClass[train_labels_oneClass!="PDAC"] = "notPDAC"
val_labels_oneClass = val_labels
val_labels_oneClass[val_labels_oneClass!="PDAC"] = "notPDAC"
#parameter search
param_res = param_search(train_counts_norm[DE_genes,],train_labels_oneClass,10,"RF")
params = list(mtry=1,replace=TRUE,nodesize=1,ntree=250) #best average parameter choices

res = classifier_suite(train_counts_norm[DE_genes,],train_labels_oneClass,meta[ids_train,],
                       val_counts_norm[DE_genes,],val_labels_oneClass,meta[ids_val,],params)
kfold_b_vs_p = res[[1]]
kfold_nc_vs_p = res[[2]]
kfold_all = res[[3]]
val_b_vs_p = res[[4]]
val_nc_vs_p = res[[5]]
val_all = res[[6]]

#### Classification when liver metastasis samples are removed ####
res = classifier_suite(train_counts_norm[DE_genes,ids_noLivMet[ids_train]],train_labels_oneClass[ids_noLivMet[ids_train]],meta[ids_noLivMet & ids_train,],
                       val_counts_norm[DE_genes,ids_noLivMet[ids_val]],val_labels_oneClass[ids_noLivMet[ids_val]],meta[ids_noLivMet & ids_val,],params)
kfold_noLiv_b_vs_p = res[[1]]
kfold_noLiv_nc_vs_p = res[[2]]
kfold_noLiv_all = res[[3]]
val_noLiv_b_vs_p = res[[4]]
val_noLiv_nc_vs_p = res[[5]]
val_noLiv_all = res[[6]]


#### Classification comparison with CA 19-9 biomarker ####
#Cross validation and external validation, using just CA 19-9 values
res = classifier_suite(rbind(1,meta$CA19_9)[,ids_train & ids_ca19],train_labels_oneClass[ids_ca19[ids_train]],meta[ids_ca19 & ids_train,],
                       rbind(1,meta$CA19_9)[,ids_val & ids_ca19],val_labels_oneClass[ids_ca19[ids_val]],meta[ids_ca19 & ids_val,],params)
kfold_ca19_b_vs_p = res[[1]]
kfold_ca19_nc_vs_p = res[[2]]
kfold_ca19_all = res[[3]]
val_ca19_b_vs_p = res[[4]]
val_ca19_nc_vs_p = res[[5]]
val_ca19_all = res[[6]]

#Cross validation and external validation, using DE genes for samples that have CA 19-9 measurements
res = classifier_suite(train_counts_norm[DE_genes,ids_ca19[ids_train]],train_labels_oneClass[ids_ca19[ids_train]],meta[ids_ca19 & ids_train,],
                       val_counts_norm[DE_genes,ids_ca19[ids_val]],val_labels_oneClass[ids_ca19[ids_val]],meta[ids_ca19 & ids_val,],params)
kfold_DE_b_vs_p = res[[1]]
kfold_DE_nc_vs_p = res[[2]]
kfold_DE_all = res[[3]]
val_DE_b_vs_p = res[[4]]
val_DE_nc_vs_p = res[[5]]
val_DE_all = res[[6]]

#Cross validation and external validation, using DE genes and CA 19-9 values
DE_ca19 = rbind(rna_counts_norm[DE_genes,],CA19_9=meta$CA19_9)
res = classifier_suite(DE_ca19[,ids_ca19 & ids_train],train_labels_oneClass[ids_ca19[ids_train]],meta[ids_ca19 & ids_train,],
                       DE_ca19[,ids_ca19 & ids_val],val_labels_oneClass[ids_ca19[ids_val]],meta[ids_ca19 & ids_val,],params)
kfold_DEca19_b_vs_p = res[[1]]
kfold_DEca19_nc_vs_p = res[[2]]
kfold_DEca19_all = res[[3]]
val_DEca19_b_vs_p = res[[4]]
val_DEca19_nc_vs_p = res[[5]]
val_DEca19_all = res[[6]]


#### Survival Analysis ####
quans = c(.5,.5,.5)
survival_data = read.csv("PDAC_survival.csv")
survival_data$Source[ survival_data$Source == "BCC_2018"] = "BCC_2019"
scores_all = c(kfold_all[[2]][,"PDAC"],val_all[[2]][,"PDAC"])
scores_pvnc = c(kfold_nc_vs_p[[2]][,"PDAC"],val_nc_vs_p[[2]][,"PDAC"])
scores_pvb = c(kfold_b_vs_p[[2]][,"PDAC"],val_b_vs_p[[2]][,"PDAC"])
score_list = list(scores_pvb,scores_pvnc,scores_all)
survival_all = survival_analysis(scores_all,survival_data,"PDAC vs All",quans[3])
survival_pvnc = survival_analysis(scores_pvnc,survival_data,"PDAC vs Non-Cancer",quans[2])
survival_pvb = survival_analysis(scores_pvb,survival_data,"PDAC vs Benign",quans[1])


#### Generate figures ####
#normalization figure
normalization_figure(rna_counts,intrinsic_genes)
#deconvolution panels
deconvolution_perc_facet(meta,dc_list,list(ids_all,ids_noLivMet,ids_stage12),c("All Samples","No Liver Met","Stage 1 & 2"),dir="Figures/",tis_num=5)
atlas_heat(atlas)
#assemble into one figure
decon_figure()
#Auc plots: cross validation and external validaiton
figure_auc_facet(meta[ids_train,],list(kfold_b_vs_p,kfold_nc_vs_p,kfold_all),"AUC_train",group_colors)
figure_auc_facet(meta[ids_val,],list(val_b_vs_p,val_nc_vs_p,val_all),"AUC_val",group_colors)
figure_auc_facet(meta[ids_noLivMet & ids_train,],list(kfold_noLiv_b_vs_p,kfold_noLiv_nc_vs_p,kfold_noLiv_all),"AUC_train_noLiverMet",group_colors)
figure_auc_facet(meta[ids_noLivMet & ids_val,],list(val_noLiv_b_vs_p,val_noLiv_nc_vs_p,val_noLiv_all),"AUC_val_noLiverMet",group_colors)
#CA19-9 comparison
train_lists = list(list(kfold_ca19_b_vs_p,kfold_ca19_nc_vs_p,kfold_ca19_all),list(kfold_DE_b_vs_p,kfold_DE_nc_vs_p,kfold_DE_all),list(kfold_DEca19_b_vs_p,kfold_DEca19_nc_vs_p,kfold_DEca19_all))
test_lists = list(list(val_ca19_b_vs_p,val_ca19_nc_vs_p,val_ca19_all),list(val_DE_b_vs_p,val_DE_nc_vs_p,val_DE_all),list(val_DEca19_b_vs_p,val_DEca19_nc_vs_p,val_DEca19_all))
figure_full_ca19(meta,train_lists,test_lists,c("CA19-9","Genes","Genes+CA19-9"))
#Survival plots
figure_survival_plot(meta,score_list,survival_pvb,survival_pvnc,survival_all,quans)
#Gene expression plots
biomaker_expression(rna_counts_norm,DE_genes,meta,group_levels)
#Parameter testing plots
param_plots(param_res)
#Age and gender correlation
age_gender_correlation(score_list,meta)