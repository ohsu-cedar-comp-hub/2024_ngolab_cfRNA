library(pheatmap)
library(ggplot2)
library(magick)
library(ggpubr)
library(pROC)
library(survminer)

de_volcano_plot = function(DE,meta,dir="../Figures/"){
  de_ids = DE$FDR < .05
  colors = c("black","red")
  names(colors) = c(FALSE,TRUE)
  p=ggplot()+
          geom_point(data=data.frame(x=as.numeric(DE$LogFC),y= -log2(as.numeric(DE$FDR)),sig=de_ids),aes(x=x,y=y,color=sig),show.legend = FALSE) + 
          geom_abline(intercept = -log2(.05),slope=0,color="red")+
          geom_text(data=data.frame(x=as.numeric(DE$LogFC[de_ids]),y= -log2(as.numeric(DE$FDR[de_ids]))+.1,Gene=DE$genes[de_ids]),
                     aes(x=x,y=y,label=Gene),check_overlap = TRUE,hjust=.5,size=4)+
          xlab("Log2 FC PDAC vs Not-PDAC")+ylab("-Log2(FDR)")+scale_color_manual(values=colors)+
          theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  
  
  ggsave(paste(dir,"DE_volcano.svg",sep=""),plot=p,device="svg",height=8,width=10)
  
  df = data.frame(Gene=DE$genes,LogFC=as.numeric(DE$LogFC),FDR= as.numeric(DE$FDR),significant=de_ids)
  write.csv(df,paste0(dir,"/source_data/supp_4_DE_volcano.csv"),row.names=FALSE)
}

deconvolution_perc_facet = function(meta,dc_list,id_list,names,fname,dir="../Figures/",tis_num=5){
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
    geom_bar(position="fill",stat="identity")+xlab("")+ylab("Relative Contribution")+
    geom_hline(yintercept = .5,linetype="dashed")+scale_fill_manual(values = small_colors)+
    theme(axis.text.x = element_text(angle=35,vjust=1,hjust=1))+ggtitle("")+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    theme(text=element_text(size=7))
  
  df = data.frame(Tissue=factor(D[,1],levels=tis_names),Expression=as.numeric(D[,2]),Group=factor(D[,3],levels=small_groups),split=factor(D[,4],levels=names))
  write.csv(df,paste0(dir,"/source_data/",fname,".csv"),row.names=FALSE)
  return(p_dc)
}

atlas_heat = function(atlas,max_tissues=29,dir="../Figures/"){
  library(ComplexHeatmap)
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
  #p_enr = pheatmap(t(atlas_norm),cluster_rows = FALSE,cluster_cols = TRUE,main = "Tissue Expression (Z-score)",annotation_row = tot_enr,annotation_legend = FALSE)
  h_anno = rowAnnotation(Total=anno_barplot(tot_enr,axis_param = list(direction = "reverse"),
                         axis=FALSE,border=FALSE),show_legend=FALSE,annotation_name_gp = gpar(fontsize=6),
                         annotation_name_rot=0)
  p_enr = Heatmap(t(atlas_norm),cluster_rows = FALSE,cluster_columns = TRUE,left_annotation = h_anno,name="Z-score",
                  row_title_gp = gpar(fontsize = 5),column_title_gp = gpar(fontsize = 5),
                  row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 6)))
                  
  p_enr = grid.grabExpr(draw(p_enr,column_title="Tissue Expression",column_title_gp = gpar(fontsize = 7)))
  
  #p_tot = ggplot(data.frame(val=tot_enr[,1],gene=factor(colnames(atlas_norm),levels=rev(colnames(atlas_norm)))),aes(x=val,y=gene))+
  #  geom_bar(stat="identity",color="gray",width=1)+xlab("")+ylab("")+
  #  theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="white"))
 
  df = data.frame(t(atlas_norm))
  write.csv(df,paste0(dir,"/source_data/fig_3a_tissue_exp.csv"))
  return(p_enr)
}

decon_figure = function(p_exp,p_tdp,dir="../Figures/"){
  #compile together
  p = ggarrange(p_exp,p_tdp,nrow=1,labels=c("a","b"),widths=c(3,2),font.label = list(size = 7))
  
  #save figure
  ggsave(paste0(dir,"Fig3_biomarker_source.pdf"),plot=p,device="pdf",heigh=80,width=180,units="mm")
}

decon_figure_split = function(p_train,p_val,dir="../Figures/"){
  #compile together
  p = ggarrange(p_train,p_val,nrow=1,labels=c("a","b"),widths=c(1,1),common.legend = TRUE,legend = "right")
  
  #save figure
  ggsave(paste0(dir,"decon_split.svg"),plot=p,device="svg",heigh=6,width=10)
}


normalization_figure = function(data_raw,intrinsic_genes,dir="../Figures/"){
  p_image = image_read(paste(dir,"figure_3_diagram.png",sep=""))
  p_image = image_ggplot(p_image)
  
  data_raw = data_raw
  int_score = intrinsic_genes$Intrinsic
  bat_score = intrinsic_genes$Extrinsic
  
  data_int = cf_norm(data_raw,intrinsic_genes)
  data_tmm = tmm(cf_norm(data_raw,intrinsic_genes))
  p = ncol(data_raw)
  
  norm_levels = c("Raw Reads","cf-Normalization","cf-Normalization + TMM")
  
  Data = c()
  Data = rbind(Data,cbind(1:p,(colSums(data_raw*int_score)+1),"Intrinsic",norm_levels[1]))
  Data = rbind(Data,cbind(1:p,(colSums(data_raw*bat_score)+1),"Extrinsic",norm_levels[1]))
  Data = rbind(Data,cbind(1:p,(colSums(data_int*int_score)+1),"Intrinsic",norm_levels[2]))
  Data = rbind(Data,cbind(1:p,(colSums(data_int*bat_score)+1),"Extrinsic",norm_levels[2]))
  Data = rbind(Data,cbind(1:p,(colSums(data_tmm*int_score)+1),"Intrinsic",norm_levels[3]))
  Data = rbind(Data,cbind(1:p,(colSums(data_tmm*bat_score)+1),"Extrinsic",norm_levels[3]))
  
  colors = c("olivedrab3","goldenrod1")
  names(colors) = c("Intrinsic","Extrinsic")
  
  p_norms = ggplot(data.frame(Contribution=factor(Data[,3],levels=names(colors)),vals=as.numeric(Data[,2]),Sample=as.numeric(Data[,1]),Norm=factor(Data[,4],levels=norm_levels)),aes(color=Contribution,x=Sample,y=vals,fill=Contribution))+
    geom_bar(position="stack", stat="identity")+xlab("Sample")+ylab("Total Counts")+facet_grid(cols=vars(Norm))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
    theme(text=element_text(size=7))
  
  p = ggarrange(p_image,p_norms,ncol=1,nrow=2,heights = c(1,.75),labels=c("a","b"),font.label = list(size = 7))
  ggsave(paste0(dir,"Fig2_normalization.pdf"),plot=p,device="pdf",height=135,width = 180, units="mm")
  
  df = data.frame(Contribution=factor(Data[,3],levels=names(colors)),value=as.numeric(Data[,2]),Sample=as.numeric(Data[,1]),Norm=factor(Data[,4],levels=norm_levels))
  write.csv(df,paste0(dir,"/source_data/fig_2b_normalization.csv"),row.names=FALSE)
}

#Create faceted auc plot for a set of 3 classifiers: P vs B, P vs NC, P vs All
#Only pass in meta for train (kfold validation) or test (external validation) split
figure_auc_facet = function(meta,class_list,file_name,fig_name,group_colors,dir="../Figures/",age_cut = 65,main_fig=FALSE){
  if(main_fig){
    lab_size=2
    point_size=.6
  } else{
    lab_size=2.7
    point_size=1
  }
  #Load PDAC scores
  score_pvb = class_list[[1]][[2]][,"PDAC"]
  score_pvnc = class_list[[2]][[2]][,"PDAC"]
  score_all = class_list[[3]][[2]][,"PDAC"]
  
  
  #filter and order meta data
  PvB = meta$Group %in% c("PDAC","Benign")
  PvNC = ! meta$Group %in% c("Other Cancer")
  
  
  ids_all = rep(TRUE,nrow(meta))
  ids_age = list(ids_all,meta$Age.at.Collection >= age_cut, meta$Age.at.Collection < age_cut)
  names(ids_age) = c("All","66 and Older","Younger than 66")
  exp_levels = c("PDAC vs Benign","PDAC vs Non-Cancer","PDAC vs All")
  ##### Make ROC Column #####
  #Compute roc values
  roc_pvb_a = roc_splits(lapply(ids_age,function(x){x[PvB]}),factor(meta$Group[PvB]=="PDAC"),score_pvb,exp_levels[1])
  roc_pvnc_a = roc_splits(lapply(ids_age,function(x){x[PvNC]}),factor(meta$Group[PvNC]=="PDAC"),score_pvnc,exp_levels[2])
  roc_all_a = roc_splits(ids_age,factor(meta$Group=="PDAC"),score_all,exp_levels[3])
  
  color_set = c("black","seagreen","orangered","slateblue","orchid")
  names(color_set) = c("All","Male","Female",names(ids_age)[-1])
  
  #make facet plot
  roc_facet_a = rbind(roc_pvb_a[[1]],roc_pvnc_a[[1]],roc_all_a[[1]])
  
  y_lab=c(.35,.2,.05)
  auc_facet_a = cbind(rbind(roc_pvb_a[[2]],roc_pvnc_a[[2]],roc_all_a[[2]]),y_lab)
  
  p_roc_a = ggplot(data.frame(y=as.numeric(roc_facet_a[,1]),x=as.numeric(roc_facet_a[,2]),Experiment=factor(roc_facet_a[,3],levels=exp_levels),Group=roc_facet_a[,4]),
                   aes(x=x,y=y,color=Group))+
    xlab("False Positive Rate")+ylab("True Positive Rate")+ggtitle("")+
    geom_path(show.legend = FALSE)+geom_abline(slope=1,linetype="dashed")+facet_grid(rows=vars(Experiment))+
    geom_label(data=data.frame(AUC=auc_facet_a[,1],Experiment=factor(auc_facet_a[,2],levels=exp_levels),Group=auc_facet_a[,3],y=as.numeric(auc_facet_a[,4])),
               aes(x=.5,y=y,label=AUC,color=Group),size=lab_size,show.legend=FALSE)+scale_color_manual(values=color_set)+
    scale_x_continuous(breaks=c(0, .5, 1))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  if('Sex' %in% colnames(meta)){
    ids_sex = list(All = ids_all,Male= meta$Sex=="M",Female= !meta$Sex=="M")
    roc_pvb_g = roc_splits(lapply(ids_sex,function(x){x[PvB]}),factor(meta$Group[PvB]=="PDAC"),score_pvb,exp_levels[1])
    roc_pvnc_g = roc_splits(lapply(ids_sex,function(x){x[PvNC]}),factor(meta$Group[PvNC]=="PDAC"),score_pvnc,exp_levels[2])
    roc_all_g = roc_splits(ids_sex,factor(meta$Group=="PDAC"),score_all,exp_levels[3])
    roc_facet_g = rbind(roc_pvb_g[[1]],roc_pvnc_g[[1]],roc_all_g[[1]])
    auc_facet_g = cbind(rbind(roc_pvb_g[[2]],roc_pvnc_g[[2]],roc_all_g[[2]]),y_lab)
    p_roc_g = ggplot(data.frame(y=as.numeric(roc_facet_g[,1]),x=as.numeric(roc_facet_g[,2]),Experiment=factor(roc_facet_g[,3],levels=exp_levels),Group=roc_facet_g[,4]),
                     aes(x=x,y=y,color=Group))+
      xlab("False Positive Rate")+ylab("True Positive Rate")+ggtitle("")+
      geom_path(show.legend = FALSE)+geom_abline(slope=1,linetype="dashed")+facet_grid(rows=vars(Experiment))+
      geom_label(data=data.frame(AUC=auc_facet_g[,1],Experiment=factor(auc_facet_g[,2],levels=exp_levels),Group=auc_facet_g[,3],y=as.numeric(auc_facet_g[,4])),
                 aes(x=.5,y=y,label=AUC,color=Group),size=lab_size,show.legend=FALSE)+scale_color_manual(values=color_set)+
      scale_x_continuous(breaks=c(0, .5, 1))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  } else{
    p_roc_g = NULL
  }
  
  ##### Make group boxplot column #####
  gcl = group_count_labels(meta$Group,group_levels)
  meta$Group = gcl[[1]]
  group_levels = gcl[[2]]
  names(group_colors) = group_levels
  #make dataframes
  group_b_facet = c()
  group_b_facet = rbind(group_b_facet,cbind(score_pvb,meta$Group[PvB],exp_levels[1],meta$SeqID[PvB]))
  group_b_facet = rbind(group_b_facet,cbind(score_pvnc,meta$Group[PvNC],exp_levels[2],meta$SeqID[PvNC]))
  group_b_facet = rbind(group_b_facet,cbind(score_all,meta$Group,exp_levels[3],meta$SeqID))
  
  #make facet plots
  p_group_b = ggplot(data.frame(Score = as.numeric(group_b_facet[,1]),Group=factor(group_b_facet[,2],levels=group_levels),Experiment=factor(group_b_facet[,3],levels=exp_levels)),
                     aes(x=Group,y=Score,color=Group))+
    geom_boxplot(outlier.shape = NA,show.legend = F)+geom_jitter(width = .2,show.legend = F,size=point_size)+xlab("Group")+ylab("PDAC Score")+ggtitle("")+
    scale_colour_manual(values=group_colors)+xlab("")+
    facet_grid(rows=vars(Experiment))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  
  
  ##### Arrange meta plots #####
  if(main_fig){ 
    p_roc_g = p_roc_g + theme(text=element_text(size=7))
    p_roc_a = p_roc_a + theme(text=element_text(size=7))
    p_group_b = p_group_b + theme(text=element_text(size=7))
    p_b = ggarrange(p_roc_g,p_roc_a,nrow=1,widths = c(.5,.5),labels=c("b","c"),align="hv",legend="bottom",common.legend=TRUE,font.label = list(size = 7))
    p_b = ggarrange(p_group_b,p_b,nrow=1,widths=c(1,1),labels=c("a",""),font.label = list(size = 7))
    ggsave(paste(dir,fig_name,"_",file_name,".pdf",sep=""),plot=p_b,device="pdf",height=108,width=180,units="mm")
  } else{
    p_b = ggarrange(p_roc_g,p_roc_a,nrow=1,widths = c(.5,.5),labels=c("b","c"),align="hv",legend="bottom",common.legend=TRUE)
    p_b = ggarrange(p_group_b,p_b,nrow=1,widths=c(1,1),labels=c("a",""))
    ggsave(paste(dir,file_name,".svg",sep=""),plot=p_b,device="svg",height=6,width=10)
  }
  
  
  
  df = data.frame(Sample=group_b_facet[,4],Score = as.numeric(group_b_facet[,1]),Group=factor(group_b_facet[,2],levels=group_levels),Experiment=factor(group_b_facet[,3],levels=exp_levels))
  write.csv(df,paste0(dir,"/source_data/",fig_name,"a_",file_name,".csv"),row.names=FALSE)
  
  df = data.frame(TPR=as.numeric(roc_facet_a[,1]),FPR=as.numeric(roc_facet_a[,2]),Experiment=factor(roc_facet_a[,3],levels=exp_levels),Age=roc_facet_a[,4])
  write.csv(df,paste0(dir,"/source_data/",fig_name,"c_",file_name,".csv"),row.names=FALSE)
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
figure_full_ca19 = function(pdac_meta,train_lists,test_lists,g_names,fig_name,dir="../Figures/"){
  test_ids = pdac_meta$Cohort == "BCC" & !is.na(pdac_meta$CA19_9)
  train_ids = pdac_meta$Cohort == "CEDAR" & !is.na(pdac_meta$CA19_9)
  
  train_plots = figure_ca19(pdac_meta[train_ids,],train_lists[[1]],train_lists[[2]],train_lists[[3]],g_names,"CA19_9_train",fig_name,c("a","c"),main_fig=TRUE)
  test_plots = figure_ca19(pdac_meta[test_ids,],test_lists[[1]],test_lists[[2]],test_lists[[3]],g_names,"CA19_9_test",fig_name,c("b","d"),main_fig=TRUE)
  
  
  ##### Arrange meta plots #####
  p_ca19_train = train_plots[[1]]
  p_roc_train = train_plots[[3]]
  p_ca19_test = test_plots[[1]]
  p_roc_test = test_plots[[3]]
  p_ca = ggarrange(p_ca19_train,p_ca19_test,ncol=1,heights=c(1,1),labels=c("a","b"),font.label = list(size = 7))
  p_roc = ggarrange(p_roc_train,p_roc_test,nrow=1,widths = c(1,1),labels=c("c","d"),align="hv",legend="bottom",common.legend=TRUE,font.label = list(size = 7))
  p_b = ggarrange(p_ca,p_roc,nrow=1,widths=c(1,1))
  #save
  #ggsave(paste(dir,fig_name,"_CA19_9_comparison.svg",sep=""),plot=p_b,device="svg",height=6,width=9)
  ggsave(paste(dir,fig_name,"_CA19_9_comparison.pdf",sep=""),plot=p_b,device="pdf",height=120,width=180,units="mm")
  
  return(c(train_plots[[4]],test_plots[[4]]))
}

#Create faceted auc plot for 2 groups of 3 classifiers: P vs B, P vs NC, P vs All
#Only pass in meta for train or test split
#pdac_figure_ca19(pdac_meta[test_ids,],list(class_ca_PvB,class_ca_PvNC,class_ca_all),list(class_gene_PvB,class_gene_PvNC,class_gene_all),list(class_both_PvB,class_both_PvNC,class_both_all),c("CA19-9","Genes","Genes+CA19-9"),"CA19_9_val")
#
figure_ca19 = function(meta,class_list1,class_list2,class_list3,g_names,name,fig_name,panel_names,dir="../Figures/",main_fig=FALSE){
  if(main_fig){
    lab_size = 2
  } else{
    lab_size = 2.7
  }
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
  
  #test significance
  sig_b = roc.test(roc_b_pvb1,roc_b_pvb2,method="delong")
  sig_nc = roc.test(roc_b_pvnc1,roc_b_pvnc2,method="delong")
  sig_all = roc.test(roc_b_all1,roc_b_all2,method="delong")
  
  
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
  
  if(sig_b$p.value < .05){ auc_b2[1] = paste0(auc_b2[1],"*") }
  if(sig_nc$p.value < .05){ auc_b2[2] = paste0(auc_b2[2],"*") }
  if(sig_all$p.value < .05){ auc_b2[3] = paste0(auc_b2[3],"*") }
  
  print(c(sig_b$p.value,sig_nc$p.value,sig_all$p.value))
  
  p_roc_b = ggplot(data.frame(y=as.numeric(roc_b_facet[,1]),x=as.numeric(roc_b_facet[,2]),Experiment=factor(roc_b_facet[,3],levels=exp_levels),Method=roc_b_facet[,4]),aes(x=x,y=y,color=Method))+
    xlab("False Positive Rate")+ylab("True Positive Rate")+ggtitle("")+
    geom_path(aes(color=Method))+geom_abline(slope=1,linetype="dashed")+facet_grid(rows=vars(Experiment))+
    geom_label(data=data.frame(AUC=auc_b1,Experiment=factor(exp_levels,levels=exp_levels)),aes(x=.5,y=.05,label=AUC),size=lab_size,show.legend = FALSE,color=set_colors[1])+
    geom_label(data=data.frame(AUC=auc_b2,Experiment=factor(exp_levels,levels=exp_levels)),aes(x=.5,y=.2,label=AUC),size=lab_size,show.legend = FALSE,color=set_colors[2])+
    geom_label(data=data.frame(AUC=auc_b3,Experiment=factor(exp_levels,levels=exp_levels)),aes(x=.5,y=.35,label=AUC),size=lab_size,show.legend = FALSE,color=set_colors[3])+
    scale_x_continuous(breaks=c(0, .5, 1))+theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
    scale_color_manual(values=set_colors)
  
  if(main_fig){
    p_roc_b = p_roc_b + theme(text=element_text(size=7))
  }
  
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
  p_ca19 = figure_ppv_npv(meta,log(meta$CA19_9),log(37),"Log CA19-9",paste0(fig_name,panel_names[1]),show_counts = TRUE,main_fig=main_fig)
  NPV = sum(meta$Group!="PDAC" & meta$CA19_9 < 37) / sum(meta$CA19_9 < 37)
  
  df = data.frame(TPR=as.numeric(roc_b_facet[,1]),FPR=as.numeric(roc_b_facet[,2]),Experiment=factor(roc_b_facet[,3],levels=exp_levels),Method=roc_b_facet[,4])
  write.csv(df,paste0(dir,"/source_data/",fig_name,panel_names[2],"_",name,".csv"),row.names=FALSE)
  
  return(list(p_ca19,p_group_b,p_roc_b,NPV))
}

group_count_labels = function(labels,groups){
  for(i in 1:length(groups)){
    ids = (labels == groups[i])
    groups[i] = paste0(groups[i],"\n(",sum(ids),")")
    labels[ids] = groups[i]
  }
  return(list(labels,groups))
}

figure_pdac_ppv_npv = function(meta,pdac_scores,ids_train,ids_val,name,score_name,fig_name,dir="../Figures/",npvs=NULL,cutoff=NULL){
  if(is.null(npvs)){ npvs = c(.95,.95) }
  if(is.null(cutoff)){
    train_NPV_cut = NPV_cutoff_value(pdac_scores[ids_train],meta$Group[ids_train]=="PDAC",NPV_target = npvs[1])
    val_NPV_cut = NPV_cutoff_value(pdac_scores[ids_val],meta$Group[ids_val]=="PDAC",NPV_target = npvs[2])
  } else{
    train_NPV_cut = cutoff
    val_NPV_cut = cutoff
  }
  fig_pdac_ppv_train = figure_ppv_npv(meta[ids_train,],pdac_scores[ids_train],train_NPV_cut,score_name,paste0(fig_name,"c"))
  fig_pdac_ppv_val = figure_ppv_npv(meta[ids_val,],pdac_scores[ids_val],val_NPV_cut,score_name,paste0(fig_name,"d"))
  
  p_ppv = ggarrange(fig_pdac_ppv_train,fig_pdac_ppv_val,ncol=1,heights=c(1,1),labels=c("c","d"))
  
  ## Without other cancer
  o_ids = meta$Group=="Other Cancer"
  if(is.null(cutoff)){
    train_NPV_cut = NPV_cutoff_value(pdac_scores[ids_train & !o_ids],meta$Group[ids_train & !o_ids]=="PDAC")
    val_NPV_cut = NPV_cutoff_value(pdac_scores[ids_val & !o_ids],meta$Group[ids_val & !o_ids]=="PDAC")
  } else{
    train_NPV_cut = cutoff
    val_NPV_cut = cutoff
  }
  
  fig_pdac_ppv_train = figure_ppv_npv(meta[ids_train & !o_ids,],pdac_scores[ids_train & !o_ids],train_NPV_cut,score_name,paste0(fig_name,"a"))
  fig_pdac_ppv_val = figure_ppv_npv(meta[ids_val & !o_ids,],pdac_scores[ids_val & !o_ids],val_NPV_cut,score_name,paste0(fig_name,"b"))
  
  p2_ppv = ggarrange(fig_pdac_ppv_train,fig_pdac_ppv_val,ncol=1,heights=c(1,1),labels=c("a","b"))
  
  p_full = ggarrange(p2_ppv,p_ppv,ncol=2)
  #save
  ggsave(paste(dir,name,"_PPV_NPV.svg",sep=""),plot=p_full,device="svg",height=6,width=9)
  
}

figure_ppv_npv = function(meta,scores,cutoff,score_name,fig_name,show_counts=TRUE,dir="../Figures/",main_fig=FALSE){
  if(main_fig){
    lab_size = 2
  } else{
    lab_size = 2.7
  }
  PPV = sum(meta$Group=="PDAC" & scores >= cutoff) / sum(scores >= cutoff)
  NPV = sum(meta$Group!="PDAC" & scores < cutoff) / sum(scores < cutoff)
  
  group_colors = c("forestgreen","blue","gold","purple","orange","red")
  group_levels = c("Benign","Pancreatitis","IPMN","PDAC","ICT","Other Cancer")
  if(show_counts){
    gcl = group_count_labels(meta$Group,group_levels)
    meta$Group = gcl[[1]]
    group_levels = gcl[[2]]
  }
  names(group_colors) = group_levels
  
   p_ppv = ggplot(data.frame(Group=factor(meta$Group,levels=group_levels),Level=scores),aes(x=Group,y=Level))+
    geom_boxplot(aes(color=Group),outlier.shape = NA,show.legend = F)+geom_hline(yintercept=cutoff,linetype="dashed",color="black")+
    geom_jitter(aes(color=Group),width=.2,show.legend=FALSE)+
    xlab("")+ylab(score_name)+ggtitle("")+
    scale_color_manual(values=group_colors)+
    geom_label(aes(label=paste0("PPV = ",round(PPV,3)),x=group_levels[2],y=max(scores)),size=lab_size)+
    geom_label(aes(label=paste0("NPV = ",round(NPV,3)),x=group_levels[3],y=max(scores)),size=lab_size)+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))
  if(main_fig){
    p_ppv = p_ppv + theme(text=element_text(size=7))
  }
  
  df = data.frame(Sample=meta$SeqID,Group=factor(meta$Group,levels=group_levels),Level=scores)
  write.csv(df,paste0(dir,"/source_data/",fig_name,"_PPV_NPV.csv"),row.names=FALSE)
   
  return(p_ppv) 
}

stage_single_boxplot = function(meta,source,score,fig_name,dir="../Figures/"){
  meta$Stage[meta$Group != "PDAC"] = ""
  meta$Stage[meta$Group == "Benign"] = "Benign"
  meta$Stage[meta$Stage == "1"] = "1-2"
  meta$Stage[meta$Stage == "2"] = "1-2"
  meta$Stage[meta$Stage == "3"] = "3-4"
  meta$Stage[meta$Stage == "4"] = "3-4"
  
  #filter and order meta data
  ids = match(names(score),meta$SeqID)
  meta = meta[ids,]
  
  if(!is.null(source)){
    all = meta$Cohort==source
  } else{
    all = rep(TRUE,nrow(meta))
  }
  
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
  
  stage_c_facet[stage_c_facet[,2]!="Benign",2] = paste("Stage",stage_c_facet[stage_c_facet[,2]!="Benign",2])
  df = data.frame(PDAC_Score = as.numeric(stage_c_facet[,1]),Stage=stage_c_facet[,2])
  write.csv(df,paste0(dir,"/source_data/",fig_name,"_boxplot.csv"),row.names=FALSE)
  return(p_stage_c)
}

figure_survival_plot_coxph = function(pdac_metadata,score_list,surv_data_pvb,surv_data_pvnc,surv_data_all,quans,dir="../Figures/"){
  #create survival plots 
  #plots that keep the distributions so it is easy to compare
  p_surv_pvb_c =  survival_single_plot_coxph(  surv_data_pvb[[1]][[2]],  surv_data_pvb[[1]][[3]],  surv_data_pvb[[1]][[4]],quans[1],NULL)
  p_surv_pvnc_c = survival_single_plot_coxph( surv_data_pvnc[[1]][[2]], surv_data_pvnc[[1]][[3]], surv_data_pvnc[[1]][[4]],quans[2],NULL)
  p_surv_all_c =  survival_single_plot_coxph(  surv_data_all[[1]][[2]],  surv_data_all[[1]][[3]],  surv_data_all[[1]][[4]],quans[3],"supp_21b")
  p_surv_pvb_b =  survival_single_plot_coxph(  surv_data_pvb[[2]][[2]],  surv_data_pvb[[2]][[3]],  surv_data_pvb[[2]][[4]],quans[1],NULL)
  p_surv_pvnc_b = survival_single_plot_coxph( surv_data_pvnc[[2]][[2]], surv_data_pvnc[[2]][[3]], surv_data_pvnc[[2]][[4]],quans[2],NULL)
  p_surv_all_b =  survival_single_plot_coxph(  surv_data_all[[2]][[2]],  surv_data_all[[2]][[3]],  surv_data_all[[2]][[4]],quans[3],"supp_21d")
  p_surv_pvb =  survival_single_plot_coxph(  surv_data_pvb[[3]][[2]],  surv_data_pvb[[3]][[3]],  surv_data_pvb[[3]][[4]],quans[1],NULL)
  p_surv_pvnc = survival_single_plot_coxph( surv_data_pvnc[[3]][[2]], surv_data_pvnc[[3]][[3]], surv_data_pvnc[[3]][[4]],quans[2],NULL)
  p_surv_all =  survival_single_plot_coxph(  surv_data_all[[3]][[2]],  surv_data_all[[3]][[3]],  surv_data_all[[3]][[4]],quans[3],"supp_21f")
  
  #compute coefficient table
  coef_table = c()
  coef_table = cbind(p_surv_pvb_c[[2]][,c(1,2,6)],cohort="CEDAR",experiment="PDAC vs Benign")
  coef_table = rbind(coef_table,cbind(p_surv_pvnc_c[[2]][,c(1,2,6)],cohort="CEDAR",experiment="PDAC vs Non-Cancer"))
  coef_table = rbind(coef_table,cbind(p_surv_all_c[[2]][,c(1,2,6)],cohort="CEDAR",experiment="PDAC vs All"))
  coef_table = rbind(coef_table,cbind(p_surv_pvb_b[[2]][,c(1,2,6)],cohort="BCC",experiment="PDAC vs Benign"))
  coef_table = rbind(coef_table,cbind(p_surv_pvnc_b[[2]][,c(1,2,6)],cohort="BCC",experiment="PDAC vs Non-Cancer"))
  coef_table = rbind(coef_table,cbind(p_surv_all_b[[2]][,c(1,2,6)],cohort="BCC",experiment="PDAC vs All"))
  coef_table = rbind(coef_table,cbind(p_surv_pvb[[2]][,c(1,2,6)],cohort="All",experiment="PDAC vs Benign"))
  coef_table = rbind(coef_table,cbind(p_surv_pvnc[[2]][,c(1,2,6)],cohort="All",experiment="PDAC vs Non-Cancer"))
  coef_table = rbind(coef_table,cbind(p_surv_all[[2]][,c(1,2,6)],cohort="All",experiment="PDAC vs All"))
  
  
  write.csv(coef_table,paste0(dir,"/source_data/supp_table6_survival_coefficients.csv"),row.names=FALSE)
  
  #Boxplots for stages
  #p_box_pvb_c = stage_single_boxplot(pdac_metadata,"CEDAR",score_list[[1]])
  #p_box_pvnc_c = stage_single_boxplot(pdac_metadata,"CEDAR",score_list[[2]])
  p_box_all_c = stage_single_boxplot(pdac_metadata,"CEDAR",score_list[[3]],"supp_21a")
  #p_box_pvb_b = stage_single_boxplot(pdac_metadata,"BCC",score_list[[1]])
  #p_box_pvnc_b = stage_single_boxplot(pdac_metadata,"BCC",score_list[[2]])
  p_box_all_b = stage_single_boxplot(pdac_metadata,"BCC",score_list[[3]],"supp_21c")
  #p_box_pvb = stage_single_boxplot(pdac_metadata,NULL,score_list[[1]])
  #p_box_pvnc = stage_single_boxplot(pdac_metadata,NULL,score_list[[2]])
  p_box_all = stage_single_boxplot(pdac_metadata,NULL,score_list[[3]],"supp_21e")
  
  
  #Arrange columns for final figures
  #wd = c(1.,2.)
  wd = c(1.,1.)
  ### Other plots that can be pulled from ###
  #p_pvb_c = ggarrange(p_box_pvb_c,p_surv_pvb_c[[1]], nrow=1,widths = wd,labels=c("a","b"))
  #p_pvb_c = annotate_figure(p_pvb_c,top=text_grob("CEDAR Cross Validation, PDAC vs Benign",size=16))
  #p_pvb_b = ggarrange(p_box_pvb_b,p_surv_pvb_b[[1]], nrow=1,widths = wd,labels=c("c","d"))
  #p_pvb_b = annotate_figure(p_pvb_b,top=text_grob("BCC Validation, PDAC vs Benign",size=16))
  #p_pvnc_c = ggarrange(p_box_pvnc_c,p_surv_pvnc_c[[1]], nrow=1, labels=c("a","b"),widths = wd)
  #p_pvnc_c = annotate_figure(p_pvnc_c,top=text_grob("CEDAR Cross Validation, PDAC vs Non-Cancer",size=16))
  #p_all_c = ggarrange(p_box_all_c,p_surv_all_c[[1]], nrow=1, labels=c("c","d"),widths = wd)
  #p_all_c = annotate_figure(p_all_c,top=text_grob("CEDAR Cross Validation, PDAC vs All",size=16))
  #p_pvnc_b = ggarrange(p_box_pvnc_b,p_surv_pvnc_b[[1]], nrow=1, labels=c("e","f"),widths = wd)
  #p_pvnc_b = annotate_figure(p_pvnc_b,top=text_grob("BCC Validation, PDAC vs Non-Cancer",size=16))
  #p_all_b = ggarrange(p_box_all_b,p_surv_all_b[[1]], nrow=1, labels=c("g","h"),widths = wd)
  #p_all_b = annotate_figure(p_all_b,top=text_grob("BCC Validation, PDAC vs All",size=16))
  #p_pvb = ggarrange(p_box_pvb,p_surv_pvb[[1]], nrow=1,widths = wd,labels=c("a","b"))
  #p_pvb = annotate_figure(p_pvb,top=text_grob("PDAC vs Benign",size=16))
  #p_pvnc = ggarrange(p_box_pvnc,p_surv_pvnc[[1]], nrow=1, labels=c("c","d"),widths = wd)
  #p_pvnc = annotate_figure(p_pvnc,top=text_grob("PDAC vs Non-Cancer",size=16))
  
  ### Supplemental figure panels ###
  p_all_c = ggarrange(p_box_all_c,p_surv_all_c[[1]], nrow=1, labels=c("a","b"),widths = wd)
  p_all_c = annotate_figure(p_all_c,top=text_grob("CEDAR Cross Validation",size=16))
  p_all_b = ggarrange(p_box_all_b,p_surv_all_b[[1]], nrow=1, labels=c("c","d"),widths = wd)
  p_all_b = annotate_figure(p_all_b,top=text_grob("BCC Validation",size=16))
  p_all = ggarrange(p_box_all,p_surv_all[[1]], nrow=1, labels=c("e","f"),widths = wd)
  p_all = annotate_figure(p_all,top=text_grob("All Cohorts",size=16))
  
  p_sup = ggarrange(p_all_c,p_all_b,p_all,ncol=1)
  
  ggsave(paste(dir,"supp_survival_coxph.svg",sep=""),plot=p_sup,device="svg",height=18,width=14.)
  
}

survival_single_plot_coxph = function(pval_in, surv_obj, surv_data_in, quan,fig_name,dir="../Figures/"){
  col_groups = c("blue","red")
  #now we want to calculate the cox proportional hazards for this model as well
  quan2 <- quan
  proportional_hazards_discrete <-surv_fit(surv_obj~(Discriminant.Score>quan2),data=surv_data_in)
  pval <- pval_in[[1]]
  print(pval)
  proportional_hazards_continuous <- coxph(surv_obj~Discriminant.Score, data=surv_data_in)
  proportional_hazards_continuous_age_sex <-  coxph(surv_obj~ Discriminant.Score+sex+age+PDAC.Stage , data = surv_data_in)
  names(col_groups) = c(paste0("Score < ",100*quan,"th percentile"),paste0("Score > ",100*quan,"th percentile"))
  p = ggsurvplot(proportional_hazards_discrete, conf.int = TRUE, legend.labs=names(col_groups),
                 ggtheme = theme_bw(), pval=FALSE, legend= c(0.3, 0.15))
  p$plot <- p$plot + annotate("label", label=pval, x=30.0,y=0.9)
  
  
  #now we make the table and combine the plot with the table
  coxph_table<- as.data.frame(summary(proportional_hazards_continuous_age_sex)$coefficients) 
  coxph_table = cbind(var=row.names(coxph_table),coxph_table)
  row.names(coxph_table) = NULL
  
  df = cbind(surv_data_in[,c(1,2,4,5,10)],Score_high=surv_data_in$Discriminant.Score>quan2)
  if(!is.null(fig_name)){
    write.csv(df,paste0(dir,"/source_data/",fig_name,"_survival.csv"),row.names=FALSE)
  }
  return(list(p$plot, coxph_table))
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
  
  df = param
  colnames(df) = c("Split_Vars","Replacement","Leaf_node_size","Number_trees","Avg_AUC")
  write.csv(df,paste0(dir,"/source_data/supp_9_parameter_search.csv"),row.names=FALSE)
}

biomaker_expression = function(data,genes,meta,group_levels,rows=5,cols=8,dir="../Figures/"){
  bcc = meta$Cohort == "BCC"
  cedar = meta$Cohort == "CEDAR"
  
  b_markers = list()
  c_markers = list()
  df = c()
  for(i in 1:length(genes)){
    b_markers[[i]] = ggplot(data.frame(val = data[genes[i],bcc],Group=factor(meta$Group[bcc],levels=group_levels)),aes(x=Group,y=val,color=Group))+
      geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2,show.legend = F)+xlab("")+ylab("")+ggtitle(genes[i])+
      scale_colour_manual(values=group_colors)+
      theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    
    df = rbind(df,data.frame(Gene = genes[i],Value = data[genes[i],bcc],Sample=meta$SeqID[bcc],Group=factor(meta$Group[bcc],levels=group_levels),Cohort="BCC"))
    
    c_markers[[i]] = ggplot(data.frame(val = data[genes[i],cedar],Group=factor(meta$Group[cedar],levels=group_levels)),aes(x=Group,y=val,color=Group))+
      geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2,show.legend = F)+xlab("")+ylab("")+ggtitle(genes[i])+
      scale_colour_manual(values=group_colors)+
      theme(panel.background = element_rect(fill="white"),panel.grid = element_line(color="gray90"))+
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    
    df = rbind(df,data.frame(Gene = genes[i],Value = data[genes[i],cedar],Sample=meta$SeqID[cedar],Group=factor(meta$Group[cedar],levels=group_levels),Cohort="CEDAR"))
    
  }
  
  labels = letters[1:length(genes)]
  p_c = annotate_figure(ggarrange(plotlist=c_markers,nrow=rows,ncol=cols,common.legend = TRUE),top=text_grob("CEDAR Log Expression"))
  p_b = annotate_figure(ggarrange(plotlist=b_markers,nrow=rows,ncol=cols,common.legend = TRUE),top=text_grob("BCC Log Expression"))
  ggsave(paste0(dir,"expression_bcc.svg"),plot=p_b,device="svg",height=10,width=15)
  ggsave(paste0(dir,"expression_cedar.svg"),plot=p_c,device="svg",height=10,width=15)
  
  write.csv(df,paste0(dir,"/source_data/supp_5_6_biomarker_expression.csv"),row.names=FALSE)
  
}

smote_comparison = function(X_sm,Y_sm,Xt,dir="../Figures/"){
  library(umap)
  library(ggsci)
  X_sm = X_sm[,Y_sm == "PDAC"]
  ### Identify synthetic samples and PDAC labels
  synth_ids = which(!colnames(X_sm) %in% colnames(Xt))
  Y = rep("Real PDAC CEDAR",ncol(X_sm))
  Y[synth_ids] = "SMOTE PDAC"
  ### Visualize PCA and umap landscape of synthetic and real PDAC samples
  um = umap(t(X_sm))
  pca = prcomp(X_sm)
  
  p_pca = dim_plot(pca$rotation[,1],pca$rotation[,2],Y,"Source","PCA of PDAC Samples",lab="PC")+
          scale_color_startrek()
  p_um = dim_plot(um$layout[,1],um$layout[,2],Y,"Source","UMAP of PDAC Samples",lab="UMAP")+
          scale_color_startrek()
  
  p = ggarrange(p_pca,p_um,nrow=1,labels = c("a","b"),common.legend = TRUE,legend="right")
  ggsave(paste0(dir,"smote_comp.svg"),plot=p,device="svg",height=6,width=12)
  
  df = data.frame(Type=Y,PCA_1=pca$rotation[,1],PCA_2=pca$rotation[,2],UMAP_1=um$layout[,1],UMAP_2=um$layout[,2])
  write.csv(df,paste0(dir,"/source_data/supp_10_smote_compare.csv"),row.names=FALSE)
}

dim_plot = function(x,y,meta,meta_name,title,lab="UMAP"){
  xrng = c(min(x),max(x))
  xadd = (xrng[2]-xrng[1])*.05
  yrng = c(min(y),max(y))
  yadd = (yrng[2]-yrng[1])*.05
  
  p = ggplot(data.frame(x=as.numeric(x),y=as.numeric(y),color=meta),
             aes(x=x,y=y,color=color))+geom_point()+
    theme(panel.background = element_rect(fill="white"),panel.grid = element_blank())+
    theme(axis.line = element_blank(),
          axis.title = element_text(vjust = 0, hjust = 0),
          axis.text=element_blank(),
          axis.ticks=element_blank())+
    annotate(geom="segment",x=(xrng[1]-xadd),xend=(xrng[1]+4*xadd),y=(yrng[1]-yadd),yend=(yrng[1]-yadd),
             arrow=grid::arrow(length=unit(.2,"cm")))+
    annotate(geom="segment",x=(xrng[1]-xadd),xend=(xrng[1]-xadd),y=(yrng[1]-yadd),yend=(yrng[1]+4*yadd),
             arrow=grid::arrow(length=unit(.2,"cm")))+
    ggtitle(title)+xlab(paste0(lab,"_1"))+ylab(paste0(lab,"_2"))+
    guides(color=guide_legend(title=meta_name))
  return(p)
}
