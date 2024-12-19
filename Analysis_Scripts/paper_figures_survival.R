require(survminer)
library(survminer)
require(ggpmisc)
library(ggpmisc)
figure_survival_plot_coxph = function(pdac_metadata,score_list,surv_data_pvb,surv_data_pvnc,surv_data_all,quans,dir="../Figures/"){
  #create survival plots 
  #plots that keep the distributions so it is easy to compare
  p_surv_pvb_c =  survival_single_plot_coxph(  surv_data_pvb[[1]][[2]],  surv_data_pvb[[1]][[3]],  surv_data_pvb[[1]][[4]],quans[1])
  p_surv_pvnc_c = survival_single_plot_coxph( surv_data_pvnc[[1]][[2]], surv_data_pvnc[[1]][[3]], surv_data_pvnc[[1]][[4]],quans[2])
  p_surv_all_c =  survival_single_plot_coxph(  surv_data_all[[1]][[2]],  surv_data_all[[1]][[3]],  surv_data_all[[1]][[4]],quans[3])
  p_surv_pvb_b =  survival_single_plot_coxph(  surv_data_pvb[[2]][[2]],  surv_data_pvb[[2]][[3]],  surv_data_pvb[[2]][[4]],quans[1])
  p_surv_pvnc_b = survival_single_plot_coxph( surv_data_pvnc[[2]][[2]], surv_data_pvnc[[2]][[3]], surv_data_pvnc[[2]][[4]],quans[2])
  p_surv_all_b =  survival_single_plot_coxph(  surv_data_all[[2]][[2]],  surv_data_all[[2]][[3]],  surv_data_all[[2]][[4]],quans[3])
 
  #Boxplots for stages
  p_box_pvb_c = stage_single_boxplot(pdac_metadata,"CEDAR",score_list[[1]])
  p_box_pvnc_c = stage_single_boxplot(pdac_metadata,"CEDAR",score_list[[2]])
  p_box_all_c = stage_single_boxplot(pdac_metadata,"CEDAR",score_list[[3]])
  p_box_pvb_b = stage_single_boxplot(pdac_metadata,"BCC",score_list[[1]])
  p_box_pvnc_b = stage_single_boxplot(pdac_metadata,"BCC",score_list[[2]])
  p_box_all_b = stage_single_boxplot(pdac_metadata,"BCC",score_list[[3]])
  
  
  #Arrange columns for final figures
  #wd = c(1.,2.)
  wd = c(1.,1.,1.5)
  p_pvb_c = ggarrange(p_box_pvb_c,p_surv_pvb_c[[1]], p_surv_pvb_c[[2]], nrow=1,widths = wd,labels=c("a","b","c"))
  p_pvb_c = annotate_figure(p_pvb_c,top=text_grob("CEDAR Cross Validation, PDAC vs Benign",size=16))
  p_pvb_b = ggarrange(p_box_pvb_b,p_surv_pvb_b[[1]], p_surv_pvb_b[[2]], nrow=1,widths = wd,labels=c("c","d","e"))
  p_pvb_b = annotate_figure(p_pvb_b,top=text_grob("BCC Validation, PDAC vs Benign",size=16))
  
  p_pvnc_c = ggarrange(p_box_pvnc_c,p_surv_pvnc_c[[1]], p_surv_pvnc_c[[2]],nrow=1, labels=c("a","b", "c"),widths = wd)
  p_pvnc_c = annotate_figure(p_pvnc_c,top=text_grob("CEDAR Cross Validation, PDAC vs Non-Cancer",size=16))
  p_all_c = ggarrange(p_box_all_c,p_surv_all_c[[1]], p_surv_all_c[[2]],nrow=1, labels=c("c","d", "e"),widths = wd)
  p_all_c = annotate_figure(p_all_c,top=text_grob("CEDAR Cross Validation, PDAC vs All",size=16))
  p_pvnc_b = ggarrange(p_box_pvnc_b,p_surv_pvnc_b[[1]], p_surv_pvnc_b[[2]],nrow=1, labels=c("f","g", "h"),widths = wd)
  p_pvnc_b = annotate_figure(p_pvnc_b,top=text_grob("BCC Validation, PDAC vs Non-Cancer",size=16))
  p_all_b = ggarrange(p_box_all_b,p_surv_all_b[[1]], p_surv_all_b[[2]],nrow=1, labels=c("i","j", "k"),widths = wd)
  p_all_b = annotate_figure(p_all_b,top=text_grob("BCC Validation, PDAC vs All",size=16))
  
  p_fig = ggarrange(p_pvb_c,p_pvb_b,nrow=2)
  p_sup = ggarrange(p_pvnc_c,p_all_c,p_pvnc_b,p_all_b,ncol=1)
  
  ggsave(paste(dir,"survival_PvB_coxph.svg",sep=""),plot=p_fig,device="svg",height=10,width=14.)
  ggsave(paste(dir,"supp_survival_coxph.svg",sep=""),plot=p_sup,device="svg",height=18,width=14.)
  
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


survival_single_plot_coxph = function(pval_in, surv_obj, surv_data_in, quan){
  col_groups = c("blue","red")
 #now we want to calculate the cox proportional hazards for this model as well
  quan2 <- quan
  proportional_hazards_discrete <-surv_fit(surv_obj~(Discriminant.Score>quan2),data=surv_data_in)
  pval <- pval_in[[1]]
  print(pval)
  proportional_hazards_continuous <- coxph(surv_obj~Discriminant.Score, data=surv_data_in)
  proportional_hazards_continuous_age_sex <-  coxph(surv_obj~ Discriminant.Score+sex+age , data = surv_data_in)
  names(col_groups) = c(paste0("Score < ",100*quan,"th percentile"),paste0("Score > ",100*quan,"th percentile"))
  p = ggsurvplot(proportional_hazards_discrete, conf.int = TRUE, legend.labs=names(col_groups),
                 ggtheme = theme_bw(), pval=FALSE, legend= c(0.3, 0.15))
  p$plot <- p$plot + annotate("label", label=pval, x=30.0,y=0.9)

  #now we make the table and combine the plot with the table
  coefficients_table<- as.data.frame(summary(proportional_hazards_continuous_age_sex)$coefficients) 
  selected_statistics_table <- coefficients_table[,c(1,3,5)] %>%
    apply(., 1, function(x) round(x, digits=3)) %>% t() %>% as.data.frame()
  selected_statistics_table$name<-c("Score", "Sex Male", "Age")
  colnames(selected_statistics_table) <- c("coef", "se(coef)", "Pr(>|z|)", "Name")
  selected_statistics_table = selected_statistics_table[, c("Name","coef", "se(coef)", "Pr(>|z|)")]
  coxph_table_grob <- tableGrob(selected_statistics_table,
                                theme=ttheme_minimal())
  coxph_table_out<- ggplot() +theme_void()+
    ylim(c(0,1)) +
    xlim(c(0,1)) +
    annotate(geom="table",
    y = .8 ,x=0.0,
    size = 18/.pt,
    label = list(data.frame(selected_statistics_table))
    ) +
    annotate(geom="text", 
             label="Cox Proportional Hazards Model",
             x=0.45,y=0.9,
             size = 20/.pt)
   
  #  geom_table(selected_statistics_table)
  #+ annotation_custom(coxph_table_grob)+
  #  theme_bw() +labs(title=" Cox Regression Model")

  combined_output <- grid.arrange(p$plot, coxph_table_grob, nrow=1)
  return(list(p$plot, coxph_table_out, combined_output,proportional_hazards_discrete, proportional_hazards_continuous_age_sex))
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
  
  all = meta$Cohort==source
  
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

