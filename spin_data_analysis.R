umap_plot = function(um,meta,meta_name,title){
  xrng = c(min(um$layout[,1]),max(um$layout[,1]))
  xadd = (xrng[2]-xrng[1])*.05
  yrng = c(min(um$layout[,2]),max(um$layout[,2]))
  yadd = (yrng[2]-yrng[1])*.05
  
  p = ggplot(data.frame(x=as.numeric(um$layout[,1]),y=as.numeric(um$layout[,2]),color=meta),
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
      ggtitle(title)+xlab("UMAP_1")+ylab("UMAP_2")+
      guides(color=guide_legend(title=meta_name))
  
  return(p)
}

source("normalization.R")
library(stringr)
### Load Spin Data ###
spin_data = read.table("../../cfRNA/Data/Spin/spin_conditions_gen27.txt")
colnames(spin_data) = str_replace(colnames(spin_data),"p","P")
spin_meta = read.table("../../cfRNA/Data/Spin/spin_conditions_MD_Age_Sex.txt")
spin_data = spin_data[,spin_meta$Coded.Name]
spin_data_cpm = cpm(spin_data)

### Umap plots with metadata ###
library(umap)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
set.seed(2)
um = umap(t(spin_data_cpm))

spin_meta$Spin.Type[spin_meta$Spin.Type=="SS1"] = "Single 1,000 RCF"
spin_meta$Spin.Type[spin_meta$Spin.Type=="SS2"] = "Single 1,600 RCF"
spin_meta$Spin.Type[spin_meta$Spin.Type=="DS11"] = "Double 1,000/2,500 RCF"
spin_meta$Spin.Type[spin_meta$Spin.Type=="DS12"] = "Double 1,000/15,000 RCF"

p_spin = umap_plot(um,spin_meta$Spin.Type,"Spin","")+scale_color_startrek()
p_time = umap_plot(um,spin_meta$Time,"Time","")+scale_color_jco()
p_tube = umap_plot(um,spin_meta$EDTA.Streck,"Tube","")+scale_color_igv()
p_donor = umap_plot(um,spin_meta$Sample.Number,"Donor","")+scale_color_d3()
p_sex = umap_plot(um,spin_meta$Sex,"Sex","")+scale_color_simpsons()

p_ref = ggarrange(p_tube,p_time,p_spin,p_donor,nrow=2,ncol=2,labels=c("a","b","c","d"))
ggsave("../Figures/reference_umap.svg",plot=p_ref,device="svg",heigh=8,width=10)

df = data.frame(Sample=spin_meta$Coded.Name,UMAP_1=um$layout[,1],UMAP_2=um$layout[,2],Spin_Type=spin_meta$Spin.Type,
                Time=spin_meta$Time,Tube=spin_meta$EDTA.Streck,Donor=spin_meta$Sample.Number,Sex=spin_meta$Sex)
write.csv(df,"../Figures/source_data/supp_1_reference_umap.csv",row.names=FALSE)


### Fit NMF model ###
source("intrinsic_factors.R")
int_fact_res = intrinsic_factors(spin_data_cpm,2,5)
fact = int_fact_res[[1]]
fact = fact/rowSums(fact)
### Identify intrinsic vs batch factors ###
gene_platelet = c("PF4","PPBP","MT-ND4","MT-CO3","TAGLN2","TM2D2","CA2","IFRD1")
gene_intrin = c("ALB","APOB","C3","CP","EEF1A1","RPL26","RPS21","CTSB")
p_sum = colSums(fact[gene_platelet,])
i_sum = colSums(fact[gene_intrin,])
fact = fact[,order(p_sum)]
colnames(fact) = c("Intrinsic","Batch")

### Save intrinsic gene list ###
fact = fact[order(fact[,2]),]
count_avg = rowMeans(spin_data_cpm)
fact = cbind(fact,spin_cpm=count_avg[row.names(fact)])
write.csv(fact,"spin_data_intrinsic.csv")

### Adjust factors based on data ###
fact = as.data.frame(fact)
fact = fact[!is.na(fact$Intrinsic),]
res_cedar = adjust_factors(cpm(rna_counts[,ids_train]),as.matrix(fact[,c(1,2)]),5)
fact_cedar = res_cedar[[1]]
colnames(fact_cedar) = colnames(fact)[1:2]
fact_cedar = fact_cedar/rowSums(fact_cedar)
write.csv(fact_cedar,"updated_intrinsic.csv")

