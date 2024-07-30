library(Glimma)
library(limma)
library(DESeq2)

project_id = snakemake@params[['project_id']]

rds = snakemake@input[['rds']]
rds = snakemake@input[['list_rds']]
rds = snakemake@input[['ERCC_rds']]

cat(sprintf(c('RDS object: ',rds,'\n')))


out_path = file.path(getwd(),'results','diffexp')
dir.create(out_path)

rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))
glMDSPlot(rds, top = 1000, groups=groups.df,path=out_path,html=paste(project_id,'mds_plot',sep='.'),launch=FALSE)

##Now with genelist normalized counts
rds = readRDS(list_rds)
groups.df = as.data.frame(colData(rds))
glMDSPlot(rds, top = 1000, groups=groups.df,path=out_path,html=paste(project_id,'mds_plot_listnorm',sep='.'),launch=FALSE)

##Now with ERCC Normalized counts
rds = readRDS(ERCC_rds)
groups.df = as.data.frame(colData(rds))
glMDSPlot(rds, top = 1000, groups=groups.df,path=out_path,html=paste(project_id,'mds_plot_ERCCnorm',sep='.'),launch=FALSE)
