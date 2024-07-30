library(DESeq2)
library(ggplot2)
library(reshape2)
library(data.table)
library(plyr)
library(RColorBrewer)


rld <- snakemake@input[['rld']]
list_rld <- snakemake@input[['list_rld']]
ERCC_rld <- snakemake@input[['ERCC_rld']]
cat(sprintf(c('rld: ', rld, '\n')))


condition <- snakemake@params[['linear_model']]
cat(sprintf(c('condition: ', condition, '\n')))

project_id <- snakemake@params[['project_id']]

density_plot <- snakemake@output[['density']]
list_density_plot <- snakemake@output[['list_density']]
ERCC_density_plot <- snakemake@output[['ERCC_density']]
cat(sprintf(c('Density plot : ', density_plot, '\n')))


colors <- snakemake@params['colors']
discrete <- snakemake@params['discrete']

# function to grab the ggplot2 colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

rld = readRDS(rld)
normed_values = assay(rld)
normed_t = t(normed_values)
meta = colData(rld)


if(colors[[1]] !='NA' & discrete[[1]] =='NA'){
    if(brewer.pal.info[colors[[1]],]$maxcolors >= length(levels(meta[[condition]]))) {
        pal <- brewer.pal(length(levels(meta[[condition]])),name=colors[[1]])
    } 
} else if(discrete[[1]] != 'NA' & length(discrete)==length(levels(meta[[condition]]))){
        pal <- unlist(discrete)
} else {
        pal <- gg_color_hue(length(levels(meta[[condition]])))
}

joined_counts = cbind(meta[condition],normed_t)

x = as.data.table(joined_counts)
mm <- melt(x,id=condition)

mu <- ddply(mm, condition, summarise, grp.mean=mean(value))
pdf(density_plot)
p<-ggplot(mm, aes_string(x='value', color=condition)) +
  geom_density()+
  geom_vline(data=mu, aes_string(xintercept='grp.mean', color=condition),
             linetype="dashed") + xlab('regularized log expression') + 
  scale_color_manual(values=pal) +
  ggtitle(eval(project_id)) + theme(plot.title = element_text(hjust = 0.5))
p
dev.off()

# Gene List normalized 
#list normed
list_rld = readRDS(list_rld)
list_normed_values = assay(list_rld)
list_normed_t = t(list_normed_values)
meta = colData(list_rld)

if(colors[[1]] !='NA' & discrete[[1]] =='NA'){
    if(brewer.pal.info[colors[[1]],]$maxcolors >= length(levels(meta[[condition]]))) {
        pal <- brewer.pal(length(levels(meta[[condition]])),name=colors[[1]])
    } 
} else if(discrete[[1]] != 'NA' & length(discrete)==length(levels(meta[[condition]]))){
        pal <- unlist(discrete)
} else {
        pal <- gg_color_hue(length(levels(meta[[condition]])))
}

joined_counts = cbind(meta[condition],list_normed_t)

x = as.data.table(joined_counts)
mm <- melt(x,id=condition)

mu <- ddply(mm, condition, summarise, grp.mean=mean(value))
pdf(list_density_plot)
p<-ggplot(mm, aes_string(x='value', color=condition)) +
  geom_density()+
  geom_vline(data=mu, aes_string(xintercept='grp.mean', color=condition),
             linetype="dashed") + xlab('regularized log expression') + 
  scale_color_manual(values=pal) +
  ggtitle(eval(project_id)) + theme(plot.title = element_text(hjust = 0.5))
p
dev.off()

# ERCC normalized 
#ERCC Normed
ERCC_rld = readRDS(ERCC_rld)
ERCC_normed_values = assay(ERCC_rld)
ERCC_normed_t = t(ERCC_normed_values)
meta = colData(ERCC_rld)

if(colors[[1]] !='NA' & discrete[[1]] =='NA'){
    if(brewer.pal.info[colors[[1]],]$maxcolors >= length(levels(meta[[condition]]))) {
        pal <- brewer.pal(length(levels(meta[[condition]])),name=colors[[1]])
    } 
} else if(discrete[[1]] != 'NA' & length(discrete)==length(levels(meta[[condition]]))){
        pal <- unlist(discrete)
} else {
        pal <- gg_color_hue(length(levels(meta[[condition]])))
}

joined_counts = cbind(meta[condition],ERCC_normed_t)

x = as.data.table(joined_counts)
mm <- melt(x,id=condition)

mu <- ddply(mm, condition, summarise, grp.mean=mean(value))
pdf(ERCC_density_plot)
p<-ggplot(mm, aes_string(x='value', color=condition)) +
  geom_density()+
  geom_vline(data=mu, aes_string(xintercept='grp.mean', color=condition),
             linetype="dashed") + xlab('regularized log expression') + 
  scale_color_manual(values=pal) +
  ggtitle(eval(project_id)) + theme(plot.title = element_text(hjust = 0.5))
p
dev.off()
