#### Duplicate analysis without using intrinsic normalization ####
source("normalization.R")
source("DE_wilcox.R")
source("paper_figures.R")
source("classification_functions.R")

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
ids_train = meta$Source=="CEDAR_2020"
ids_val = meta$Source == "BCC_2019"

#### Perform normalization ####
rna_counts_tpm = tpm(rna_counts)

### Separate into training and validation sets ####
train_counts_tpm = rna_counts_tpm[,ids_train]
val_counts_tpm = rna_counts_tpm[,ids_val]
train_labels = meta$Group[ids_train]
val_labels = meta$Group[ids_val]

#remove low variance genes, computed on training set
V = apply(train_counts_tpm,1,var)
train_counts_tpm = train_counts_tpm[V > 1,]
val_counts_tpm = val_counts_tpm[V > 1,]

#### Feature selection and filtering ####
DE = DE_wilcox(train_counts_tpm,train_labels=="PDAC")
DE_tpm_genes = DE[DE[,"FDR"]<.05,1]

#### Classification ####
params = list(mtry=1,replace=TRUE,nodesize=1,ntree=250) #best average parameter choices
train_labels_oneClass = train_labels
train_labels_oneClass[train_labels_oneClass!="PDAC"] = "notPDAC"
val_labels_oneClass = val_labels
val_labels_oneClass[val_labels_oneClass!="PDAC"] = "notPDAC"

res = classifier_suite(train_counts_tpm[DE_tpm_genes,],train_labels_oneClass,meta[ids_train,],
                       val_counts_tpm[DE_tpm_genes,],val_labels_oneClass,meta[ids_val,],params)
kfold_tpm_b_vs_p = res[[1]]
kfold_tpm_nc_vs_p = res[[2]]
kfold_tpm_all = res[[3]]
val_tpm_b_vs_p = res[[4]]
val_tpm_nc_vs_p = res[[5]]
val_tpm_all = res[[6]]

##### AUC plots ####
figure_auc_facet(meta[ids_train,],list(kfold_tpm_b_vs_p,kfold_tpm_nc_vs_p,kfold_tpm_all),"AUC_train_tpm",group_colors)
figure_auc_facet(meta[ids_val,],list(val_tpm_b_vs_p,val_tpm_nc_vs_p,val_tpm_all),"AUC_val_tpm",group_colors)
