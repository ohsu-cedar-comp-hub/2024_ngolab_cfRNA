library("DESeq2")
library("pheatmap")
library("ggplot2")
library("ggrepel")

print('Setting parameters')

rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

rld = snakemake@input[['rld_out']]
cat(sprintf(c('RLD object: ',rld,'\n')))

list_rds = snakemake@input[['list_rds']]
cat(sprintf(c('list_RDS object: ',list_rds,'\n')))

list_rld = snakemake@input[['list_rld_out']]
cat(sprintf(c('list_RLD object: ',list_rld,'\n')))

ERCC_rds = snakemake@input[['ERCC_rds']]
cat(sprintf(c('ERCC_RDS object: ',ERCC_rds,'\n')))

ERCC_rld = snakemake@input[['ERCC_rld_out']]
cat(sprintf(c('ERCC_RLD object: ',ERCC_rld,'\n')))

Type = snakemake@params[['linear_model']]
cat(sprintf(c('Type: ',Type,'\n')))

sampleID = snakemake@params[['sample_id']]
cat(sprintf(c('Sample ID: ',sampleID,'\n')))

ma_plot = snakemake@output[['ma_plot']]
cat(sprintf(c('MA plot', ma_plot,'\n')))

out_table = snakemake@output[['table']]
cat(sprintf(c('Summary results table', out_table,'\n')))

out_gene_table = snakemake@output[['geneID_table']]

panel_ma = snakemake@output[['panel_ma']]
cat(sprintf(c('MA panel', panel_ma,'\n')))

heatmap_plot = snakemake@output[['heatmap_plot']]
cat(sprintf(c('Heatmap plot', heatmap_plot,'\n')))

var_heat = snakemake@output[['var_heat']]
cat(sprintf(c('Variance Heatmap plot', var_heat,'\n')))

pca_plot = snakemake@output[['pca_plot']]
cat(sprintf(c('PCA plot', pca_plot,'\n')))

list_ma_plot = snakemake@output[['list_ma_plot']]
cat(sprintf(c('list_MA plot', list_ma_plot,'\n')))

list_out_table = snakemake@output[['list_table']]
cat(sprintf(c('list_Summary results table', list_out_table,'\n')))

list_out_gene_table = snakemake@output[['list_geneID_table']]

list_panel_ma = snakemake@output[['list_panel_ma']]
cat(sprintf(c('list_MA panel', list_panel_ma,'\n')))

list_heatmap_plot = snakemake@output[['list_heatmap_plot']]
cat(sprintf(c('list_Heatmap plot', list_heatmap_plot,'\n')))

list_var_heat = snakemake@output[['list_var_heat']]
cat(sprintf(c('list_Variance Heatmap plot', list_var_heat,'\n')))

list_pca_plot = snakemake@output[['list_pca_plot']]
cat(sprintf(c('list_PCA plot', list_pca_plot,'\n')))

ERCC_ma_plot = snakemake@output[['ERCC_ma_plot']]
cat(sprintf(c('ERCC_MA plot', ERCC_ma_plot,'\n')))

ERCC_out_table = snakemake@output[['ERCC_table']]
cat(sprintf(c('ERCC_Summary results table', ERCC_out_table,'\n')))

ERCC_out_gene_table = snakemake@output[['ERCC_geneID_table']]

ERCC_panel_ma = snakemake@output[['ERCC_panel_ma']]
cat(sprintf(c('ERCC_MA panel', ERCC_panel_ma,'\n')))

ERCC_heatmap_plot = snakemake@output[['ERCC_heatmap_plot']]
cat(sprintf(c('ERCC_Heatmap plot', ERCC_heatmap_plot,'\n')))

ERCC_var_heat = snakemake@output[['ERCC_var_heat']]
cat(sprintf(c('Variance Heatmap plot', ERCC_var_heat,'\n')))

ERCC_pca_plot = snakemake@output[['ERCC_pca_plot']]
cat(sprintf(c('PCA plot', ERCC_pca_plot,'\n')))


labels <- snakemake@params[['pca_labels']]
cat(sprintf(c('PCA Labels: ',labels, '\n')))

cat(sprintf('Load dds DESeqTransform object'))
cat(sprintf(c('rds test', rds,'\n')))
cat(sprintf(c('rld rds test', rld,'\n')))
dds <- readRDS(rds)

cat(sprintf('Load rlog DESeqTransform object'))
rld <- readRDS(rld)

Dir <- "results/diffexp/pairwise/"

plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)

contrast <- c(Type, snakemake@params[["contrast"]])
baseline <- contrast[[3]]
target <- contrast[[2]]

upCol = "#FF9999"
downCol = "#99CCFF"
ncCol = "#CCCCCC"
adjp <- 0.01
FC <- 2

ens2geneID <- snakemake@config[['ens2geneID']]

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Pairwise PCA Plot
pdf(pca_plot)
plotPCA(rld, intgroup=labels)
dev.off()

# Pairwise PCA Plot with more than one PCA parameter
if (length(labels)>1) {
  pca_plot2 <- sub("$",paste(contrast[2],"vs",contrast[3],"twoDimensional_pca_plot.pdf", sep = "-"), Dir)
  pcaData <- plotPCA(rld, intgroup=c(labels[[1]], labels[[2]]), returnData=TRUE)
  pdf(pca_plot2, 5, 5)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]], shape=labels[[2]])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
  dev.off()
}

res <- results(dds, contrast=contrast,  independentFiltering = FALSE, cooksCutoff = Inf)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)

# MA plot - calc norm values yourself to plot with ggplot
# MA plot is log2normalized counts (averaged across all samples) vs. log2FC

# extract normalized counts to calculate values for MA plot
norm_counts <- counts(dds, normalized=TRUE)

## select up regulated genes
forPlot <- as.data.frame(res)
forPlot$log2Norm <- log2(rowMeans(norm_counts))
forPlot$Gene <- rownames(forPlot)

## Replace ensemble id's with gene id's
gene_id = read.delim(ens2geneID)

## Remove unique identifier .xx from heatmap data
forPlot$Gene <- sub("\\.[0-9]*", "", forPlot$Gene)
iv <- match(forPlot$Gene, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
## row number for all matches between these two dataframes
## Use paste to get rid of factors of this column, and just paste the value of the gene name
forPlot$Gene  <- paste(gene_id[iv, "external_gene_name"])

up <- forPlot$padj < adjp & forPlot$log2FoldChange > log2(FC)
sum(up)

## select down regulated genes
down <- forPlot$padj < adjp & forPlot$log2FoldChange < -log2(FC)
sum(down)

# Grab the top 5 up and down regulated genes to label in the volcano plot
if (sum(up)>5) {
  temp <- forPlot[up,]
  upGenesToLabel <- head(rownames(temp[order(-temp$log2FoldChange),], 5))
} else if (sum(up) %in% 1:5) {
  temp <- forPlot[up,]
  upGenesToLabel <- rownames(temp[order(-temp$log2FoldChange),])
}

if (sum(down)>5) {
  temp <- forPlot[down,]
  downGenesToLabel <- head(rownames(temp[order(temp$log2FoldChange),], 5))
} else if (sum(down) %in% 1:5) {
  temp <- forPlot[down,]
  downGenesToLabel <- rownames(temp[order(temp$log2FoldChange),])
}

forPlot$Expression <- ifelse(down, 'down',
                  ifelse(up, 'up','NS'))
forPlot$Expression <- factor(forPlot$Expression, levels=c("up","down","NS"))

# Assign colours to conditions
if (sum(up)==0 & sum(down)==0) {
  colours <- ncCol
} else if (sum(up)==0) {
  colours <- c(downCol, ncCol)
} else if (sum(down)==0) {
  colours <- c(upCol, ncCol)
} else {
  colours <- c(upCol, downCol, ncCol)
}

# Create vector for labelling the genes based on whether genes are DE or not
if (exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- c(downGenesToLabel, upGenesToLabel)
} else if (exists("downGenesToLabel") & !exists("upGenesToLabel")) {
  genesToLabel <- downGenesToLabel
} else if (!exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- upGenesToLabel
}

if (exists("genesToLabel")) {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    geom_label_repel(aes(label=ifelse(Gene %in% genesToLabel, as.character(Gene),'')),box.padding=0.1, point.padding=0.5, segment.color="gray70", show.legend=FALSE) +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
} else {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
}

# MA plot
pdf(ma_plot)
print({
  maPlot
})
dev.off()

# P-histogram
p_hist = snakemake@output[['p_hist']]
pdf(p_hist)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main='P values for genes with mean normalized count larger than 1',xlab='pvalue')
dev.off()

#panel ma plot
pdf(panel_ma)
par(mfrow=c(2,2),mar=c(2,2,1,1) +0.1)
ylim <- c(-2.5,2.5)
resGA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
mtext(resG@elementMetadata$description[[2]], outer=T, cex=.6,line=-1)
dev.off()

# Heatmap of top 50 genes
topGenes <- head(order(res$padj),50)

df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- df[,subset_cols]
}

forHeat <- assay(rld)[topGenes,]

## Remove unique identifier .xx from heatmap data
rownames(forHeat) <- sub("\\.[0-9]*", "", rownames(forHeat))
iv <- match(rownames(forHeat), gene_id$ensembl_gene_id)
head(gene_id[iv,])

## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
## row number for all matches between these two dataframes
## Use paste to get rid of factors of this column, and just paste the value of the gene name
rownames(forHeat) <- paste(gene_id[iv, "external_gene_name"])

pdf(heatmap_plot)
pheatmap(forHeat, cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes:", contrast[2], "vs", contrast[3]))
dev.off()

# Variance Heatmap
pdf(var_heat)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

## Remove unique identifier .xx from heatmap data
rownames(mat) <- sub("\\.[0-9]*", "", rownames(mat))
iv <- match(rownames(mat), gene_id$ensembl_gene_id)
head(gene_id[iv,])
rownames(mat) <- paste(gene_id[iv, "external_gene_name"])

pheatmap(mat, scale="row", annotation_col = annot,fontsize=6, main = paste("Heatmap of top 50 most variable genes:", contrast[2], "vs", contrast[3]))
dev.off()

# sort by p-value
res <- res[order(res$padj),]
write.table(as.data.frame(res), file=out_table, quote=FALSE, sep='\t')

# Export table with the geneIDs
## Remove unique identifier .xx from results table
res$GeneID <- sub("\\.[0-9]*", "", rownames(res))
iv <- match(res$GeneID, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## Paste external gene name of these ensembl IDs into this column instead
## Use paste to get rid of factors of this column, and just paste the value of the gene name
res$GeneID  <- paste(gene_id[iv, "external_gene_name"])

resGen <- cbind(res$GeneID, as.data.frame(res)[,1:6])
colnames(resGen)[1] <- "GeneID"

# Export to table
write.table(resGen, file=out_gene_table, row.names=FALSE, quote=FALSE, sep="\t")

##############################List normalized pairwise##################################################
cat(sprintf('Load list_dds DESeqTransform object'))
cat(sprintf(c('list rds test', list_rds,'\n')))
dds <- readRDS(list_rds)

cat(sprintf('Load list_rlog DESeqTransform object'))
cat(sprintf(c('list_rld_rds test', list_rld,'\n')))
rld <- readRDS(list_rld)


plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)

contrast <- c(Type, snakemake@params[["contrast"]])
baseline <- contrast[[3]]
target <- contrast[[2]]

upCol = "#FF9999"
downCol = "#99CCFF"
ncCol = "#CCCCCC"
adjp <- 0.01
FC <- 2

ens2geneID <- snakemake@config[['ens2geneID']]

parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

# Pairwise PCA Plot
pdf(list_pca_plot)
plotPCA(rld, intgroup=labels)
dev.off()

# Pairwise PCA Plot with more than one PCA parameter
if (length(labels)>1) {
  pca_plot2 <- sub("$",paste(contrast[2],"vs",contrast[3],"twoDimensional_pca_plot_list.pdf", sep = "-"), Dir)
  pcaData <- plotPCA(list_rld, intgroup=c(labels[[1]], labels[[2]]), returnData=TRUE)
  pdf(pca_plot2, 5, 5)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]], shape=labels[[2]])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
  dev.off()
}

res <- results(dds, contrast=contrast,  independentFiltering = FALSE, cooksCutoff = Inf)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)

# MA plot - calc norm values yourself to plot with ggplot
# MA plot is log2normalized counts (averaged across all samples) vs. log2FC

# extract normalized counts to calculate values for MA plot
norm_counts <- counts(dds, normalized=TRUE)

## select up regulated genes
forPlot <- as.data.frame(res)
forPlot$log2Norm <- log2(rowMeans(norm_counts))
forPlot$Gene <- rownames(forPlot)

## Replace ensemble id's with gene id's
gene_id = read.delim(ens2geneID)

## Remove unique identifier .xx from heatmap data
forPlot$Gene <- sub("\\.[0-9]*", "", forPlot$Gene)
iv <- match(forPlot$Gene, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
## row number for all matches between these two dataframes
## Use paste to get rid of factors of this column, and just paste the value of the gene name
forPlot$Gene  <- paste(gene_id[iv, "external_gene_name"])

up <- forPlot$padj < adjp & forPlot$log2FoldChange > log2(FC)
sum(up)

## select down regulated genes
down <- forPlot$padj < adjp & forPlot$log2FoldChange < -log2(FC)
sum(down)

# Grab the top 5 up and down regulated genes to label in the volcano plot
if (sum(up)>5) {
  temp <- forPlot[up,]
  upGenesToLabel <- head(rownames(temp[order(-temp$log2FoldChange),], 5))
} else if (sum(up) %in% 1:5) {
  temp <- forPlot[up,]
  upGenesToLabel <- rownames(temp[order(-temp$log2FoldChange),])
}

if (sum(down)>5) {
  temp <- forPlot[down,]
  downGenesToLabel <- head(rownames(temp[order(temp$log2FoldChange),], 5))
} else if (sum(down) %in% 1:5) {
  temp <- forPlot[down,]
  downGenesToLabel <- rownames(temp[order(temp$log2FoldChange),])
}

forPlot$Expression <- ifelse(down, 'down',
                             ifelse(up, 'up','NS'))
forPlot$Expression <- factor(forPlot$Expression, levels=c("up","down","NS"))

# Assign colours to conditions
if (sum(up)==0 & sum(down)==0) {
  colours <- ncCol
} else if (sum(up)==0) {
  colours <- c(downCol, ncCol)
} else if (sum(down)==0) {
  colours <- c(upCol, ncCol)
} else {
  colours <- c(upCol, downCol, ncCol)
}

# Create vector for labelling the genes based on whether genes are DE or not
if (exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- c(downGenesToLabel, upGenesToLabel)
} else if (exists("downGenesToLabel") & !exists("upGenesToLabel")) {
  genesToLabel <- downGenesToLabel
} else if (!exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- upGenesToLabel
}

if (exists("genesToLabel")) {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    geom_label_repel(aes(label=ifelse(Gene %in% genesToLabel, as.character(Gene),'')),box.padding=0.1, point.padding=0.5, segment.color="gray70", show.legend=FALSE) +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
} else {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
}

# MA plot
pdf(list_ma_plot)
print({
  maPlot
})
dev.off()

# P-histogram
list_p_hist = snakemake@output[['list_p_hist']]
pdf(list_p_hist)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main='P values for genes with mean normalized count larger than 1',xlab='pvalue')
dev.off()

#panel ma plot
pdf(list_panel_ma)
par(mfrow=c(2,2),mar=c(2,2,1,1) +0.1)
ylim <- c(-2.5,2.5)
resGA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
mtext(resG@elementMetadata$description[[2]], outer=T, cex=.6,line=-1)
dev.off()

# Heatmap of top 50 genes
topGenes <- head(order(res$padj),50)

df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- df[,subset_cols]
}

forHeat <- assay(rld)[topGenes,]

## Remove unique identifier .xx from heatmap data
rownames(forHeat) <- sub("\\.[0-9]*", "", rownames(forHeat))
iv <- match(rownames(forHeat), gene_id$ensembl_gene_id)
head(gene_id[iv,])

## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
## row number for all matches between these two dataframes
## Use paste to get rid of factors of this column, and just paste the value of the gene name
rownames(forHeat) <- paste(gene_id[iv, "external_gene_name"])

pdf(list_heatmap_plot)
pheatmap(forHeat, cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes:", contrast[2], "vs", contrast[3]))
dev.off()

# Variance Heatmap
pdf(list_var_heat)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

## Remove unique identifier .xx from heatmap data
rownames(mat) <- sub("\\.[0-9]*", "", rownames(mat))
iv <- match(rownames(mat), gene_id$ensembl_gene_id)
head(gene_id[iv,])
rownames(mat) <- paste(gene_id[iv, "external_gene_name"])

pheatmap(mat, scale="row", annotation_col = annot,fontsize=6, main = paste("Heatmap of top 50 most variable genes:", contrast[2], "vs", contrast[3]))
dev.off()

# sort by p-value
res <- res[order(res$padj),]
write.table(as.data.frame(res), file=list_out_table, quote=FALSE, sep='\t')

# Export table with the geneIDs
## Remove unique identifier .xx from results table
res$GeneID <- sub("\\.[0-9]*", "", rownames(res))
iv <- match(res$GeneID, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## Paste external gene name of these ensembl IDs into this column instead
## Use paste to get rid of factors of this column, and just paste the value of the gene name
res$GeneID  <- paste(gene_id[iv, "external_gene_name"])

resGen <- cbind(res$GeneID, as.data.frame(res)[,1:6])
colnames(resGen)[1] <- "GeneID"

# Export to table
write.table(resGen, file=list_out_gene_table, row.names=FALSE, quote=FALSE, sep="\t")

##############################ERCC normalized pairwise##################################################
cat(sprintf('Load ERCC_dds DESeqTransform object'))
cat(sprintf(c('ERCC_rds test', ERCC_rds,'\n')))
dds <- readRDS(ERCC_rds)

cat(sprintf('Load ERCC_rlog DESeqTransform object'))
cat(sprintf(c('ERCC_rld_rds test', ERCC_rld,'\n')))
rld <- readRDS(ERCC_rld)

plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)

contrast <- c(Type, snakemake@params[["contrast"]])
baseline <- contrast[[3]]
target <- contrast[[2]]

upCol = "#FF9999"
downCol = "#99CCFF"
ncCol = "#CCCCCC"
adjp <- 0.01
FC <- 2

ens2geneID <- snakemake@config[['ens2geneID']]

parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

# Pairwise PCA Plot
pdf(ERCC_pca_plot)
plotPCA(rld, intgroup=labels)
dev.off()

# Pairwise PCA Plot with more than one PCA parameter
if (length(labels)>1) {
  pca_plot2 <- sub("$",paste(contrast[2],"vs",contrast[3],"twoDimensional_pca_plot_ERCC.pdf", sep = "-"), Dir)
  pcaData <- plotPCA(ERCC_rld, intgroup=c(labels[[1]], labels[[2]]), returnData=TRUE)
  pdf(pca_plot2, 5, 5)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]], shape=labels[[2]])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
  dev.off()
}

res <- results(dds, contrast=contrast,  independentFiltering = FALSE, cooksCutoff = Inf)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)

# MA plot - calc norm values yourself to plot with ggplot
# MA plot is log2normalized counts (averaged across all samples) vs. log2FC

# extract normalized counts to calculate values for MA plot
norm_counts <- counts(dds, normalized=TRUE)

## select up regulated genes
forPlot <- as.data.frame(res)
forPlot$log2Norm <- log2(rowMeans(norm_counts))
forPlot$Gene <- rownames(forPlot)

## Replace ensemble id's with gene id's
gene_id = read.delim(ens2geneID)

## Remove unique identifier .xx from heatmap data
forPlot$Gene <- sub("\\.[0-9]*", "", forPlot$Gene)
iv <- match(forPlot$Gene, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
## row number for all matches between these two dataframes
## Use paste to get rid of factors of this column, and just paste the value of the gene name
forPlot$Gene  <- paste(gene_id[iv, "external_gene_name"])

up <- forPlot$padj < adjp & forPlot$log2FoldChange > log2(FC)
sum(up)

## select down regulated genes
down <- forPlot$padj < adjp & forPlot$log2FoldChange < -log2(FC)
sum(down)

# Grab the top 5 up and down regulated genes to label in the volcano plot
if (sum(up)>5) {
  temp <- forPlot[up,]
  upGenesToLabel <- head(rownames(temp[order(-temp$log2FoldChange),], 5))
} else if (sum(up) %in% 1:5) {
  temp <- forPlot[up,]
  upGenesToLabel <- rownames(temp[order(-temp$log2FoldChange),])
}

if (sum(down)>5) {
  temp <- forPlot[down,]
  downGenesToLabel <- head(rownames(temp[order(temp$log2FoldChange),], 5))
} else if (sum(down) %in% 1:5) {
  temp <- forPlot[down,]
  downGenesToLabel <- rownames(temp[order(temp$log2FoldChange),])
}

forPlot$Expression <- ifelse(down, 'down',
                             ifelse(up, 'up','NS'))
forPlot$Expression <- factor(forPlot$Expression, levels=c("up","down","NS"))

# Assign colours to conditions
if (sum(up)==0 & sum(down)==0) {
  colours <- ncCol
} else if (sum(up)==0) {
  colours <- c(downCol, ncCol)
} else if (sum(down)==0) {
  colours <- c(upCol, ncCol)
} else {
  colours <- c(upCol, downCol, ncCol)
}

# Create vector for labelling the genes based on whether genes are DE or not
if (exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- c(downGenesToLabel, upGenesToLabel)
} else if (exists("downGenesToLabel") & !exists("upGenesToLabel")) {
  genesToLabel <- downGenesToLabel
} else if (!exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- upGenesToLabel
}

if (exists("genesToLabel")) {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    geom_label_repel(aes(label=ifelse(Gene %in% genesToLabel, as.character(Gene),'')),box.padding=0.1, point.padding=0.5, segment.color="gray70", show.legend=FALSE) +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
} else {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
}

# MA plot
pdf(ERCC_ma_plot)
print({
  maPlot
})
dev.off()

# P-histogram
ERCC_p_hist = snakemake@output[['ERCC_p_hist']]
pdf(ERCC_p_hist)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main='P values for genes with mean normalized count larger than 1',xlab='pvalue')
dev.off()

#panel ma plot
pdf(ERCC_panel_ma)
par(mfrow=c(2,2),mar=c(2,2,1,1) +0.1)
ylim <- c(-2.5,2.5)
resGA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
mtext(resG@elementMetadata$description[[2]], outer=T, cex=.6,line=-1)
dev.off()

# Heatmap of top 50 genes
topGenes <- head(order(res$padj),50)

df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- df[,subset_cols]
}

forHeat <- assay(rld)[topGenes,]

## Remove unique identifier .xx from heatmap data
rownames(forHeat) <- sub("\\.[0-9]*", "", rownames(forHeat))
iv <- match(rownames(forHeat), gene_id$ensembl_gene_id)
head(gene_id[iv,])

## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
## row number for all matches between these two dataframes
## Use paste to get rid of factors of this column, and just paste the value of the gene name
rownames(forHeat) <- paste(gene_id[iv, "external_gene_name"])

pdf(ERCC_heatmap_plot)
pheatmap(forHeat, cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes:", contrast[2], "vs", contrast[3]))
dev.off()

# Variance Heatmap
pdf(ERCC_var_heat)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

## Remove unique identifier .xx from heatmap data
rownames(mat) <- sub("\\.[0-9]*", "", rownames(mat))
iv <- match(rownames(mat), gene_id$ensembl_gene_id)
head(gene_id[iv,])
rownames(mat) <- paste(gene_id[iv, "external_gene_name"])

pheatmap(mat, scale="row", annotation_col = annot,fontsize=6, main = paste("Heatmap of top 50 most variable genes:", contrast[2], "vs", contrast[3]))
dev.off()

# sort by p-value
res <- res[order(res$padj),]
write.table(as.data.frame(res), file=ERCC_out_table, quote=FALSE, sep='\t')

# Export table with the geneIDs
## Remove unique identifier .xx from results table
res$GeneID <- sub("\\.[0-9]*", "", rownames(res))
iv <- match(res$GeneID, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## Paste external gene name of these ensembl IDs into this column instead
## Use paste to get rid of factors of this column, and just paste the value of the gene name
res$GeneID  <- paste(gene_id[iv, "external_gene_name"])

resGen <- cbind(res$GeneID, as.data.frame(res)[,1:6])
colnames(resGen)[1] <- "GeneID"

# Export to table
write.table(resGen, file=ERCC_out_gene_table, row.names=FALSE, quote=FALSE, sep="\t")
