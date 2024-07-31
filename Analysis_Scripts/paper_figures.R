library(pheatmap)
library(ggplot2)
library(magick)
library(ggpubr)
library(pROC)

deconvolution_perc_facet = function(meta,dc_list,id_list,names,dir="../Figures/",tis_num=5){
  small_groups = c("Other","PDAC")
  small_colors = c("gray","purple")
  names(small_colors) = small_groups
  
  #Create dataframe of PDAC and non-PDAC tissue percentages, sorted by tissues with largest PDAC percentage
  D = c()
  tis_names = names(sort((colMeans(dc_list[[1]][meta$Group[id_list[[1]]]=="PDAC",])-colMeans(dc_list[[1]][meta$Group[id_list[[1]]]!="PDAC",]))/colMeans(dc_list[[1]][meta$Group[id_list[[1]]]!="PDAC",]),decreasing=TRUE)[1:tis_num])
  for(i in 1:length(dc_list)){
    tis_pdac = colMeans(dc_list[[i]][meta$Group[id_list[[i]]]=="PDAC",tis_names])
    tis_other = colMeans(dc_list[[i]][meta$Group[id_list[[i]]]!="PDAC",tis_names])
    tis_total = tis_pdac+tis_other
    D = rbind(D,cbind(tis_names,tis_pdac/tis_total,"PDAC",names[i]),cbind(tis_names,tis_other/tis_total,"Other",names[i]))
  }
  
  #Create faceted barplot
  p_dc = ggplot(data.frame(Tissue=factor(D[,1],levels=tis_names),Exp=as.numeric(D[,2]),Group=factor(D[,3],levels=small_groups),split=factor(D[,4],levels=names)),
                aes(x=Tissue,y=Exp,fill=Group))+facet_grid(rows=vars(split))+
    geom_bar(position="fill",stat="identity")+xlab("")+ylab("Relative Expression")+
    geom_hline(yintercept = .5,linetype="dashed")+scale_fill_manual(values = small_colors)+
    theme(axis.text.x = element_text(angle=35,size=10,vjust=1,hjust=1))+ggtitle("")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  ggsave(paste(dir,"tissue_decon_PDAC_facet.svg",sep=""),plot=p_dc,device="svg",height=6,width=4)
}

atlas_heat = function(atlas,max_tissues=29,dir="../Figures/"){
  #Create z-norm, top expression tissue list
  atlas_norm = atlas - rowMeans(atlas) 
  var = sqrt(apply(atlas_norm,1,var))
  var[var==0] = 1
  atlas_norm = atlas_norm / var
  
  atlas_norm[atlas_norm > -min(atlas_norm)] = -min(atlas_norm)
  
  rel_exp = colSums(atlas_norm)
  sids = sort(rel_exp,decreasing=TRUE,index.return=TRUE)$ix[1:max_tissues]
  atlas = atlas_norm
  atlas_norm = atlas_norm[,sids]
  
  #heatmap of tissue enrichment for each gene
  tot_enr = data.frame("Total" = colSums(atlas_norm))
  p_enr = pheatmap(t(atlas_norm),cluster_rows = FALSE,cluster_cols = TRUE,main = "Tissue Expression (Z-score)",annotation_row = tot_enr,annotation_legend = FALSE)
  
  p_tot = ggplot(data.frame(val=tot_enr[,1],gene=factor(colnames(atlas_norm),levels=rev(colnames(atlas_norm)))),aes(x=val,y=gene))+
    geom_bar(stat="identity",color="gray",width=1)+xlab("")+ylab("")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="white"))
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  #theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
  
  #arrange plots
  #p = ggarrange(p_enr,p_q,nrow=2,heights = c(1,.25))
  ggsave(paste(dir,"tissue_heatmap_exp.png",sep=""),plot=p_enr,device="png",height=6,width=6)
  ggsave(paste(dir,"tissue_exp_bar.png",sep=""),plot=p_tot,device="png",height=6,width=6)
}

decon_figure = function(dir="../Figures/"){
  #Load figures
  p_exp = image_read(paste(dir,"tissue_heatmap_exp.png",sep=""))
  p_exp = image_ggplot(p_exp)
  p_tdp = image_read(paste(dir,"tissue_decon_PDAC_facet.png",sep=""))
  p_tdp = image_ggplot(p_tdp)
  
  
  #compile together
  p = ggarrange(p_exp,p_tdp,nrow=1,labels=c("a","b"),widths=c(3,2))
  
  #save figure
  ggsave(paste0(dir,"biomarker_source.svg"),plot=p,device="svg",heigh=6,width=10)
}


normalization_figure = function(data_raw,int_score,dir="../Figures/"){
  p_image = image_read(paste(dir,"figure_3_diagram.png",sep=""))
  p_image = image_ggplot(p_image)
  
  #Select top k genes by score
  k=1000
  int_genes = int_score[1:100,1]
  bat_genes = int_score[nrow(int_score):(nrow(int_score)-k),1]
  
  data_int = t(t(data_raw)/colSums(data_raw[int_genes,])) * mean(colSums(data_raw[int_genes,]))
  data_tmm = tmm(intrinsic_norm(data_raw,int_genes))
  p = ncol(data_raw)
  
  norm_levels = c("Raw Reads","Intrinsic Norm","Intrinsic Norm + TMM")
  
  Data = c()
  Data = rbind(Data,cbind(1:p,log(colSums(data_raw[int_genes,])+1),"Intrinsic",norm_levels[1]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_raw[bat_genes,])+1),"Batch",norm_levels[1]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_raw)),"All",norm_levels[1]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_int[int_genes,])+1),"Intrinsic",norm_levels[2]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_int[bat_genes,])+1),"Batch",norm_levels[2]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_int)),"All",norm_levels[2]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_tmm[int_genes,])+1),"Intrinsic",norm_levels[3]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_tmm[bat_genes,])+1),"Batch",norm_levels[3]))
  Data = rbind(Data,cbind(1:p,log(colSums(data_tmm)),"All",norm_levels[3]))
  
  colors = c("skyblue","olivedrab3","goldenrod1")
  names(colors) = c("All","Intrinsic","Batch")
  
  p_norms = ggplot(data.frame(Genes=factor(Data[,3],levels=names(colors)),vals=as.numeric(Data[,2]),Sample=as.numeric(Data[,1]),Norm=factor(Data[,4],levels=norm_levels)),aes(color=Genes,x=Sample,y=vals))+
    geom_line()+xlab("Sample")+ylab("Log Total Expression")+facet_grid(cols=vars(Norm))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+scale_color_manual(values=colors)
  
  p = ggarrange(p_image,p_norms,ncol=1,nrow=2,heights = c(1,.5),labels=c("a","b"))
  ggsave(paste0(dir,"normalization.svg"),plot=p,device="svg",height=6,width=8)
}

#Create faceted auc plot for a set of 3 classifiers: P vs B, P vs NC, P vs All
#Only pass in meta for train (kfold validation) or test (external validation) split
figure_auc_facet = function(meta,class_list,file_name,group_colors,dir="../Figures/",age_cut = 68){
  #Load PDAC scores
  score_pvb = class_list[[1]][[2]][,"PDAC"]
  score_pvnc = class_list[[2]][[2]][,"PDAC"]
  score_all = class_list[[3]][[2]][,"PDAC"]
  
  #filter and order meta data
  PvB = meta$Group %in% c("PDAC","Benign")
  PvNC = ! meta$Group %in% c("Other Cancer")
  
  ids_all = rep(TRUE,nrow(meta))
  ids_gender = list(All = ids_all,Male= meta$Gender.Flag=="M",Female= !meta$Gender.Flag=="M")
  ids_age = list(ids_all,meta$Age.at.Collection >= age_cut, meta$Age.at.Collection < age_cut)
  names(ids_age) = c("All","68 and Older","Yonger than 68")
  exp_levels = c("PDAC vs Benign","PDAC vs Non-Cancer","PDAC vs All")
  ##### Make ROC Column #####
  #Compute roc values
  roc_pvb_g = roc_splits(lapply(ids_gender,function(x){x[PvB]}),factor(meta$Group[PvB]=="PDAC"),score_pvb,exp_levels[1])
  roc_pvb_a = roc_splits(lapply(ids_age,function(x){x[PvB]}),factor(meta$Group[PvB]=="PDAC"),score_pvb,exp_levels[1])
  roc_pvnc_g = roc_splits(lapply(ids_gender,function(x){x[PvNC]}),factor(meta$Group[PvNC]=="PDAC"),score_pvnc,exp_levels[2])
  roc_pvnc_a = roc_splits(lapply(ids_age,function(x){x[PvNC]}),factor(meta$Group[PvNC]=="PDAC"),score_pvnc,exp_levels[2])
  roc_all_g = roc_splits(ids_gender,factor(meta$Group=="PDAC"),score_all,exp_levels[3])
  roc_all_a = roc_splits(ids_age,factor(meta$Group=="PDAC"),score_all,exp_levels[3])
  
  
  color_set = c("black","seagreen","orangered","slateblue","orchid")
  names(color_set) = c(names(ids_gender),names(ids_age)[-1])
  
  #make facet plot
  roc_facet_g = rbind(roc_pvb_g[[1]],roc_pvnc_g[[1]],roc_all_g[[1]])
  roc_facet_a = rbind(roc_pvb_a[[1]],roc_pvnc_a[[1]],roc_all_a[[1]])
  y_lab=c(.35,.2,.05)
  auc_facet_g = cbind(rbind(roc_pvb_g[[2]],roc_pvnc_g[[2]],roc_all_g[[2]]),y_lab)
  auc_facet_a = cbind(rbind(roc_pvb_a[[2]],roc_pvnc_a[[2]],roc_all_a[[2]]),y_lab)
  
  p_roc_g = ggplot(data.frame(y=as.numeric(roc_facet_g[,1]),x=as.numeric(roc_facet_g[,2]),Experiment=factor(roc_facet_g[,3],levels=exp_levels),Group=roc_facet_g[,4]),
    aes(x=x,y=y,color=Group))+
    xlab("False Positive Rate")+ylab("True Positive Rate")+ggtitle("")+
    geom_path(show.legend = FALSE)+geom_abline(slope=1,linetype="dashed")+facet_grid(rows=vars(Experiment))+
    geom_label(data=data.frame(AUC=auc_facet_g[,1],Experiment=factor(auc_facet_g[,2],levels=exp_levels),Group=auc_facet_g[,3],y=as.numeric(auc_facet_g[,4])),
    aes(x=.5,y=y,label=AUC,color=Group),size=2.7,show.legend=FALSE)+scale_color_manual(values=color_set)+
    scale_x_continuous(breaks=c(0, .5, 1))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  p_roc_a = ggplot(data.frame(y=as.numeric(roc_facet_a[,1]),x=as.numeric(roc_facet_a[,2]),Experiment=factor(roc_facet_a[,3],levels=exp_levels),Group=roc_facet_a[,4]),
                   aes(x=x,y=y,color=Group))+
    xlab("False Positive Rate")+ylab("True Positive Rate")+ggtitle("")+
    geom_path(show.legend = FALSE)+geom_abline(slope=1,linetype="dashed")+facet_grid(rows=vars(Experiment))+
    geom_label(data=data.frame(AUC=auc_facet_a[,1],Experiment=factor(auc_facet_a[,2],levels=exp_levels),Group=auc_facet_a[,3],y=as.numeric(auc_facet_a[,4])),
               aes(x=.5,y=y,label=AUC,color=Group),size=2.7,show.legend=FALSE)+scale_color_manual(values=color_set)+
    scale_x_continuous(breaks=c(0, .5, 1))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  
  ##### Make group boxplot column #####
  #make dataframes
  group_b_facet = c()
  group_b_facet = rbind(group_b_facet,cbind(score_pvb,meta$Group[PvB],exp_levels[1]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvnc,meta$Group[PvNC],exp_levels[2]))
  group_b_facet = rbind(group_b_facet,cbind(score_all,meta$Group,exp_levels[3]))
  
  #make facet plots
  p_group_b = ggplot(data.frame(Score = as.numeric(group_b_facet[,1]),Group=factor(group_b_facet[,2],levels=group_levels),Experiment=factor(group_b_facet[,3],levels=exp_levels)),
                     aes(x=Group,y=Score,color=Group))+
    geom_boxplot(outlier.shape = NA,show.legend = F)+geom_jitter(width = .2,show.legend = F)+xlab("Group")+ylab("PDAC Score")+ggtitle("")+
    scale_colour_manual(values=group_colors)+xlab("")+
    facet_grid(rows=vars(Experiment))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  
  ##### Arrange meta plots #####
  p_b = ggarrange(p_roc_g,p_roc_a,nrow=1,widths = c(.5,.5),labels=c("b","c"),align="hv",legend="bottom",common.legend=TRUE)
  p_b = ggarrange(p_group_b,p_b,nrow=1,widths=c(1,1),labels=c("a",""))
  #save
  ggsave(paste(dir,file_name,".svg",sep=""),plot=p_b,device="svg",height=6,width=10)
}

#Compute roc values for a list of different sample sets
#Compile into a table for use in ggplot
roc_splits = function(ids_list,labels,values,exp){
  roc_df = c()
  auc_df = c()
  for(i in 1:length(ids_list)){
    rc = roc(factor(labels[ids_list[[i]]]),values[ids_list[[i]]])
    roc_df = rbind(roc_df,cbind(rc$sensitivities,1-rc$specificities,exp,names(ids_list)[i]))
    auc_df = rbind(auc_df,c(paste0(names(ids_list)[i],": ",round(rc$auc,3)),exp,names(ids_list)[i]))
  }
  return(list(roc_df,auc_df))
}

#Combine training and testing plots for CA19-9 analysis
#train_lists and test_lists are lists of 3 lists, each sublist has the PvB,PvNC,all breakdown
#train_lists = list(list(kfold_ca_PvB,kfold_ca_PvNC,kfold_ca_all),list(kfold_gene_PvB,kfold_gene_PvNC,kfold_gene_all),list(kfold_both_PvB,kfold_both_PvNC,kfold_both_all))
#test_lists = list(list(class_ca_PvB,class_ca_PvNC,class_ca_all),list(class_gene_PvB,class_gene_PvNC,class_gene_all),list(class_both_PvB,class_both_PvNC,class_both_all))
#pdac_figure_full_ca19(pdac_meta,train_lists,test_lists,c("CA19-9","Genes","Genes+CA19-9"),"CA19")
figure_full_ca19 = function(pdac_meta,train_lists,test_lists,g_names,dir="../Figures/"){
  test_ids = pdac_meta$Source %in% c("BCC_2019","BCC_2018") & !is.na(pdac_meta$CA19_9)
  train_ids = pdac_meta$Source == "CEDAR_2020" & !is.na(pdac_meta$CA19_9)
  
  train_plots = figure_ca19(pdac_meta[train_ids,],train_lists[[1]],train_lists[[2]],train_lists[[3]],g_names,"CA19_9_train")
  test_plots = figure_ca19(pdac_meta[test_ids,],test_lists[[1]],test_lists[[2]],test_lists[[3]],g_names,"CA19_9_test")
  
  
  ##### Arrange meta plots #####
  p_ca19_train = train_plots[[1]]
  p_roc_train = train_plots[[3]]
  p_ca19_test = test_plots[[1]]
  p_roc_test = test_plots[[3]]
  p_ca = ggarrange(p_ca19_train,p_ca19_test,ncol=1,heights=c(1,1),labels=c("a","b"))
  p_roc = ggarrange(p_roc_train,p_roc_test,nrow=1,widths = c(1,1),labels=c("c","d"),align="hv",legend="bottom",common.legend=TRUE)
  p_b = ggarrange(p_ca,p_roc,nrow=1)
  #save
  ggsave(paste(dir,"CA19_9_comparison.svg",sep=""),plot=p_b,device="svg",height=6,width=9)
  
}

#Create faceted auc plot for 2 groups of 3 classifiers: P vs B, P vs NC, P vs All
#Only pass in meta for train or test split
#pdac_figure_ca19(pdac_meta[test_ids,],list(class_ca_PvB,class_ca_PvNC,class_ca_all),list(class_gene_PvB,class_gene_PvNC,class_gene_all),list(class_both_PvB,class_both_PvNC,class_both_all),c("CA19-9","Genes","Genes+CA19-9"),"CA19_9_val")
#
figure_ca19 = function(meta,class_list1,class_list2,class_list3,g_names,name){
  #Load scores for training and BCC sets
  score_pvb1 = class_list1[[1]][[2]][,"PDAC"]
  score_pvnc1 = class_list1[[2]][[2]][,"PDAC"]
  score_all1 = class_list1[[3]][[2]][,"PDAC"]
  score_pvb2 = class_list2[[1]][[2]][,"PDAC"]
  score_pvnc2 = class_list2[[2]][[2]][,"PDAC"]
  score_all2 = class_list2[[3]][[2]][,"PDAC"]
  score_pvb3 = class_list3[[1]][[2]][,"PDAC"]
  score_pvnc3 = class_list3[[2]][[2]][,"PDAC"]
  score_all3 = class_list3[[3]][[2]][,"PDAC"]
  
  
  #filter and order meta data
  PvB = meta$Group %in% c("PDAC","Benign")
  PvNC = ! meta$Group %in% c("Other Cancer")
  
  exp_levels = c("PDAC vs Benign","PDAC vs Non-Cancer","PDAC vs All")
  ##### Make ROC Column #####
  #Compute roc values
  roc_b_all1 = roc(factor(meta$Group=="PDAC"),score_all1)
  roc_b_pvb1 = roc(factor(meta$Group[PvB]=="PDAC"),score_pvb1)
  roc_b_pvnc1 = roc(factor(meta$Group[PvNC]=="PDAC"),score_pvnc1)
  roc_b_all2 = roc(factor(meta$Group=="PDAC"),score_all2)
  roc_b_pvb2 = roc(factor(meta$Group[PvB]=="PDAC"),score_pvb2)
  roc_b_pvnc2 = roc(factor(meta$Group[PvNC]=="PDAC"),score_pvnc2)
  roc_b_all3 = roc(factor(meta$Group=="PDAC"),score_all3)
  roc_b_pvb3 = roc(factor(meta$Group[PvB]=="PDAC"),score_pvb3)
  roc_b_pvnc3 = roc(factor(meta$Group[PvNC]=="PDAC"),score_pvnc3)
  
  #make facet plot
  roc_b_facet = c()
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_pvb1$sensitivities,1-roc_b_pvb1$specificities,exp_levels[1],g_names[1]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_pvnc1$sensitivities,1-roc_b_pvnc1$specificities,exp_levels[2],g_names[1]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_all1$sensitivities,1-roc_b_all1$specificities,exp_levels[3],g_names[1]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_pvb2$sensitivities,1-roc_b_pvb2$specificities,exp_levels[1],g_names[2]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_pvnc2$sensitivities,1-roc_b_pvnc2$specificities,exp_levels[2],g_names[2]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_all2$sensitivities,1-roc_b_all2$specificities,exp_levels[3],g_names[2]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_pvb3$sensitivities,1-roc_b_pvb3$specificities,exp_levels[1],g_names[3]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_pvnc3$sensitivities,1-roc_b_pvnc3$specificities,exp_levels[2],g_names[3]))
  roc_b_facet = rbind(roc_b_facet,cbind(roc_b_all3$sensitivities,1-roc_b_all3$specificities,exp_levels[3],g_names[3]))
  
  set_colors = c("magenta","black","plum4")
  names(set_colors) = g_names
  auc_b1 = paste0(g_names[1],": ",round(c(roc_b_pvb1$auc,roc_b_pvnc1$auc,roc_b_all1$auc),3))
  auc_b2 = paste0(g_names[2],": ",round(c(roc_b_pvb2$auc,roc_b_pvnc2$auc,roc_b_all2$auc),3))
  auc_b3 = paste0(g_names[3],": ",round(c(roc_b_pvb3$auc,roc_b_pvnc3$auc,roc_b_all3$auc),3))
  
  p_roc_b = ggplot(data.frame(y=as.numeric(roc_b_facet[,1]),x=as.numeric(roc_b_facet[,2]),Experiment=factor(roc_b_facet[,3],levels=exp_levels),Method=roc_b_facet[,4]),aes(x=x,y=y,color=Method))+
    xlab("False Positive Rate")+ylab("True Positive Rate")+ggtitle("")+
    geom_path(aes(color=Method))+geom_abline(slope=1,linetype="dashed")+facet_grid(rows=vars(Experiment))+
    geom_label(data=data.frame(AUC=auc_b1,Experiment=factor(exp_levels,levels=exp_levels)),aes(x=.5,y=.05,label=AUC),size=2.7,show.legend = FALSE,color=set_colors[1])+
    geom_label(data=data.frame(AUC=auc_b2,Experiment=factor(exp_levels,levels=exp_levels)),aes(x=.5,y=.2,label=AUC),size=2.7,show.legend = FALSE,color=set_colors[2])+
    geom_label(data=data.frame(AUC=auc_b3,Experiment=factor(exp_levels,levels=exp_levels)),aes(x=.5,y=.35,label=AUC),size=2.7,show.legend = FALSE,color=set_colors[3])+
    scale_x_continuous(breaks=c(0, .5, 1))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    scale_color_manual(values=set_colors)
  
  ##### Make group boxplot column #####
  group_colors = c("forestgreen","blue","gold","purple","orange","red")
  group_levels = c("Benign","Pancreatitis","IPMN","PDAC","ICT","Other Cancer")
  names(group_colors) = group_levels
  #make dataframes
  group_b_facet = c()
  group_b_facet = rbind(group_b_facet,cbind(score_pvb1,meta$Group[PvB],exp_levels[1],g_names[1]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvnc1,meta$Group[PvNC],exp_levels[2],g_names[1]))
  group_b_facet = rbind(group_b_facet,cbind(score_all1,meta$Group,exp_levels[3],g_names[1]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvb2,meta$Group[PvB],exp_levels[1],g_names[2]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvnc2,meta$Group[PvNC],exp_levels[2],g_names[2]))
  group_b_facet = rbind(group_b_facet,cbind(score_all2,meta$Group,exp_levels[3],g_names[2]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvb3,meta$Group[PvB],exp_levels[1],g_names[3]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvnc3,meta$Group[PvNC],exp_levels[2],g_names[3]))
  group_b_facet = rbind(group_b_facet,cbind(score_all3,meta$Group,exp_levels[3],g_names[3]))
  
  #make facet plots
  p_group_b = ggplot(data.frame(Score = as.numeric(group_b_facet[,1]),Group=factor(group_b_facet[,2],levels=group_levels),
                                Experiment=factor(group_b_facet[,3],levels=exp_levels),Set=factor(group_b_facet[,4],levels=g_names)),
                     aes(x=Group,y=Score,color=Set,fill=Group))+
    geom_violin(show.legend = F,lwd=0.6)+#geom_jitter(show.legend = F,width=.1,aes(shape=Set))+
    xlab("")+ylab("PDAC Score")+ggtitle("")+
    scale_colour_manual(values=set_colors)+scale_fill_manual(values=group_colors)+
    facet_grid(rows=vars(Experiment))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  ##### Make CA19-9 plots with PPV/NPV #####
  PPV = sum(meta$Group=="PDAC" & meta$CA19_9 >= 37) / sum(meta$CA19_9 >= 37)
  NPV = sum(meta$Group!="PDAC" & meta$CA19_9 < 37) / sum(meta$CA19_9 < 37)
  p_ca19 = ggplot(data.frame(Group=factor(meta$Group,levels=group_levels),Level=log(as.numeric(meta$CA19_9))),aes(x=Group,y=Level))+
    geom_violin(aes(color=Group),show.legend = FALSE)+geom_hline(yintercept=log(37),linetype="dashed",color="black")+
    geom_point(aes(color=Group),show.legend=FALSE,position=position_jitterdodge())+
    xlab("")+ylab("Log CA19-9")+ggtitle("")+ylim(c(0,10))+
    scale_fill_manual(values=group_colors)+
    geom_label(aes(label=paste0("PPV = ",round(PPV,3)),x="Pancreatitis",y=9),size=2.7)+
    geom_label(aes(label=paste0("NPV = ",round(NPV,3)),x="IPMN",y=9),size=2.7)+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  
  return(list(p_ca19,p_group_b,p_roc_b))
}

figure_survival_plot = function(pdac_metadata,score_list,surv_data_pvb,surv_data_pvnc,surv_data_all,quans,dir="../Figures/"){
  #Create survival plots
  p_surv_pvb_c = survival_single_plot(surv_data_pvb[[1]][[1]],surv_data_pvb[[1]][[2]],quans[1])
  p_surv_pvnc_c = survival_single_plot(surv_data_pvnc[[1]][[1]],surv_data_pvnc[[1]][[2]],quans[2])
  p_surv_all_c = survival_single_plot(surv_data_all[[1]][[1]],surv_data_all[[1]][[2]],quans[3])
  p_surv_pvb_b = survival_single_plot(surv_data_pvb[[2]][[1]],surv_data_pvb[[2]][[2]],quans[1])
  p_surv_pvnc_b = survival_single_plot(surv_data_pvnc[[2]][[1]],surv_data_pvnc[[2]][[2]],quans[2])
  p_surv_all_b = survival_single_plot(surv_data_all[[2]][[1]],surv_data_all[[2]][[2]],quans[3])
  
  
  #Boxplots for stages
  p_box_pvb_c = stage_single_boxplot(pdac_metadata,"CEDAR_2020",score_list[[1]])
  p_box_pvnc_c = stage_single_boxplot(pdac_metadata,"CEDAR_2020",score_list[[2]])
  p_box_all_c = stage_single_boxplot(pdac_metadata,"CEDAR_2020",score_list[[3]])
  p_box_pvb_b = stage_single_boxplot(pdac_metadata,"BCC_2019",score_list[[1]])
  p_box_pvnc_b = stage_single_boxplot(pdac_metadata,"BCC_2019",score_list[[2]])
  p_box_all_b = stage_single_boxplot(pdac_metadata,"BCC_2019",score_list[[3]])
  
  
  #Arrange columns for final figures
  wd = c(.6,1)
  p_pvb_c = ggarrange(p_box_pvb_c,p_surv_pvb_c,nrow=1,widths = wd,labels=c("a","b"))
  p_pvb_c = annotate_figure(p_pvb_c,top=text_grob("CEDAR Cross Validation, PDAC vs Benign",size=16))
  p_pvb_b = ggarrange(p_box_pvb_b,p_surv_pvb_b,nrow=1,widths = wd,labels=c("c","d"))
  p_pvb_b = annotate_figure(p_pvb_b,top=text_grob("BCC Validation, PDAC vs Benign",size=16))
  
  p_pvnc_c = ggarrange(p_box_pvnc_c,p_surv_pvnc_c,labels=c("a","b"),widths = wd)
  p_pvnc_c = annotate_figure(p_pvnc_c,top=text_grob("CEDAR Cross Validation, PDAC vs Non-Cancer",size=16))
  p_all_c = ggarrange(p_box_all_c,p_surv_all_c,labels=c("c","d"),widths = wd)
  p_all_c = annotate_figure(p_all_c,top=text_grob("CEDAR Cross Validation, PDAC vs All",size=16))
  p_pvnc_b = ggarrange(p_box_pvnc_b,p_surv_pvnc_b,labels=c("e","f"),widths = wd)
  p_pvnc_b = annotate_figure(p_pvnc_b,top=text_grob("BCC Validation, PDAC vs Non-Cancer",size=16))
  p_all_b = ggarrange(p_box_all_b,p_surv_all_b,labels=c("g","h"),widths = wd)
  p_all_b = annotate_figure(p_all_b,top=text_grob("BCC Validation, PDAC vs All",size=16))
  
  p_fig = ggarrange(p_pvb_c,p_pvb_b,nrow=2)
  p_sup = ggarrange(p_pvnc_c,p_all_c,p_pvnc_b,p_all_b,ncol=1)
  
  ggsave(paste(dir,"survival_PvB.svg",sep=""),plot=p_fig,device="svg",height=6.5,width=6.5)
  ggsave(paste(dir,"supp_survival.svg",sep=""),plot=p_sup,device="svg",height=13,width=6.5)
  
}

survival_single_plot = function(surv_data,pvals,quan){
  col_groups = c("blue","red")
  names(col_groups) = c(paste0("Score < ",100*quan,"th percentile"),paste0("Score > ",100*quan,"th percentile"))
  p = ggplot(data.frame(surv=as.numeric(surv_data[,1]),time=as.numeric(surv_data[,2]),group = factor(surv_data[,3],levels=names(col_groups))),aes(x=time,y=surv,color=group))+
    geom_step(aes(linetype=group),size=1.5)+scale_color_manual(values=col_groups)+ylab("Survival Probability")+xlab("Survival Time (Months)")+
    ylim(c(0,1))+theme(legend.position = "bottom",legend.title = element_blank(),panel.background = element_rect(fill="white"),legend.key = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    geom_label(data=data.frame(pv = pvals[1]),aes(x=25,y=0.75,label=pv),color="black")
  
  return(p)
}

stage_single_boxplot = function(meta,source,score){
  meta$Stage[meta$Group != "PDAC"] = ""
  meta$Stage[meta$Group == "Benign"] = "Benign"
  meta$Stage[meta$Stage == "1"] = "1-2"
  meta$Stage[meta$Stage == "2"] = "1-2"
  meta$Stage[meta$Stage == "3"] = "3-4"
  meta$Stage[meta$Stage == "4"] = "3-4"
  
  #filter and order meta data
  ids = match(names(score),meta$SeqID)
  meta = meta[ids,]
  
  all = meta$Source==source
  
  ##### Make stage boxplot column #####
  stage_levels = c("Benign","1-2","3-4")
  stage_colors = c("forestgreen","mediumpurple1","magenta")
  names(stage_colors) = stage_levels
  #make dataframes
  stage_c_facet = c()
  stage_c_facet = rbind(stage_c_facet,cbind(score[all],meta$Stage[all]))
  stage_c_facet = stage_c_facet[stage_c_facet[,2] %in% stage_levels,]
  #make plots
  p_stage_c = ggplot(data.frame(Score = as.numeric(stage_c_facet[,1]),Stage=factor(stage_c_facet[,2],levels=stage_levels)),
                     aes(x=Stage,y=Score,color=Stage))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2)+xlab("PDAC Stage")+ylab("PDAC Score")+
    theme(legend.position="bottom")+theme(axis.text.x = element_text(size = 12))+xlab("PDAC Stage")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    scale_color_manual(values=stage_colors)+
    theme(legend.text = element_text(color="transparent"),legend.title = element_text(color="transparent"),legend.key = element_rect(fill="transparent"))+
    guides(color=guide_legend(override.aes = list(color="transparent",fill="transparent")))
  
  return(p_stage_c)
}


param_plots = function(param,dir="../Figures/"){
  lev_m = c("1","5","10")
  lev_t = c("250","500","1000")
  AUC = param[,5]
  
  p_mtry = ggplot(data.frame(val=factor(as.character(param[,1]),levels=lev_m),AUC=AUC),aes(x=val,y=AUC,color=val))+
              geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2)+xlab("")+ylab("Avg AUC")+
              theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
              ggtitle("Variables to Split")+theme(legend.position="none")
  
  replace = rep("Yes",nrow(param))
  replace[param[,2]==0] = "No"
  
  p_rep = ggplot(data.frame(val=replace,AUC=AUC),aes(x=val,y=AUC,color=val))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2)+xlab("")+ylab("Avg AUC")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    ggtitle("Sample with Replacement")+theme(legend.position="none")
  
  p_node = ggplot(data.frame(val=factor(as.character(param[,3]),levels=lev_m),AUC=AUC),aes(x=val,y=AUC,color=val))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2)+xlab("")+ylab("Avg AUC")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    ggtitle("Leaf Node Size")+theme(legend.position="none")
  
  p_trees = ggplot(data.frame(val=factor(as.character(param[,4]),levels=lev_t),AUC=AUC),aes(x=val,y=AUC,color=val))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2)+xlab("")+ylab("Avg AUC")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    ggtitle("Number of Trees")+theme(legend.position="none")
  
  p = ggarrange(p_mtry,p_rep,p_node,p_trees,ncol=2,nrow=2,labels=c("a","b","c","d"))
  ggsave(paste0(dir,"param_boxplots.svg"),plot=p,device="svg",height=4,width=6)
}

biomaker_expression = function(data,genes,meta,group_levels,rows=5,cols=8,dir="../Figures/"){
  bcc = meta$Source == "BCC_2018" | meta$Source == "BCC_2019"
  cedar = meta$Source == "CEDAR_2020"
  
  b_markers = list()
  c_markers = list()
  for(i in 1:length(genes)){
    b_markers[[i]] = ggplot(data.frame(val = data[genes[i],bcc],Group=factor(meta$Group[bcc],levels=group_levels)),aes(x=Group,y=val,color=Group))+
      geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2,show.legend = F)+xlab("")+ylab("")+ggtitle(genes[i])+
      scale_colour_manual(values=group_colors)+
      theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    
    c_markers[[i]] = ggplot(data.frame(val = data[genes[i],cedar],Group=factor(meta$Group[cedar],levels=group_levels)),aes(x=Group,y=val,color=Group))+
      geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2,show.legend = F)+xlab("")+ylab("")+ggtitle(genes[i])+
      scale_colour_manual(values=group_colors)+
      theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  }
  
  labels = letters[1:length(genes)]
  p_c = annotate_figure(ggarrange(plotlist=c_markers,nrow=rows,ncol=cols,common.legend = TRUE),top=text_grob("CEDAR Log Expression"))
  p_b = annotate_figure(ggarrange(plotlist=b_markers,nrow=rows,ncol=cols,common.legend = TRUE),top=text_grob("BCC Log Expression"))
  ggsave(paste0(dir,"expression_bcc.svg"),plot=p_b,device="svg",height=10,width=15)
  ggsave(paste0(dir,"expression_cedar.svg"),plot=p_c,device="svg",height=10,width=15)
}

age_gender_correlation = function(score_list,meta,dir="../Figures/"){
  meta = meta[meta$Group!="PDAC" & meta$Gender.Flag %in% c("M","F"),]
  exp_labels = c("PDAC vs Benign","PDAC vs Non-Cancer","PDAC vs All")
  df = rbind(cbind(score_list[[1]][names(score_list[[1]]) %in% meta$SeqID],meta$Age.at.Collection[meta$SeqID %in% names(score_list[[1]])],meta$Gender.Flag[meta$SeqID %in% names(score_list[[1]])],exp_labels[1]),
             cbind(score_list[[2]][names(score_list[[2]]) %in% meta$SeqID],meta$Age.at.Collection[meta$SeqID %in% names(score_list[[2]])],meta$Gender.Flag[meta$SeqID %in% names(score_list[[2]])],exp_labels[2]),
             cbind(score_list[[3]][names(score_list[[3]]) %in% meta$SeqID],meta$Age.at.Collection[meta$SeqID %in% names(score_list[[3]])],meta$Gender.Flag[meta$SeqID %in% names(score_list[[3]])],exp_labels[3]))
  p_age = ggplot(data.frame(score=as.numeric(df[,1]),age=as.numeric(df[,2]),exp=factor(df[,4],levels=exp_labels)),aes(x=age,y=score))+
                  geom_point()+xlab("Age at Collection")+ylab("PDAC Score")+
                  facet_grid(rows=vars(exp))
  p_gender = ggplot(data.frame(score=as.numeric(df[,1]),gender=df[,3],exp=factor(df[,4],levels=exp_labels)),aes(x=gender,y=score))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2)+xlab("Gender")+
    ylab("PDAC Score")+facet_grid(rows=vars(exp))
  
  p = ggarrange(p_age,p_gender,ncol=2,nrow=1,labels=c("a","b"))
  ggsave(paste0(dir,"age_gender.svg"),plot=p,device="svg",height=6,width=6)
}