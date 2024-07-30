library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(reshape2)
library(gridExtra)
library(grid)
library(stringr)
library(pheatmap)
library(readxl)
library(cowplot)
# output files
exon_fraction <- snakemake@output[["exon_fraction"]]
input_reads <- snakemake@output[["input_reads"]]
unique_reads <- snakemake@output[["unique_reads"]]
lib_size <- snakemake@output[["lib_size"]]
md_stats <- snakemake@output[["md_stats"]]
batch_pca <- snakemake@output[["batch_pca"]]
sub_counts <- snakemake@output[["sub_counts"]]
subset_cts <- snakemake@output[["subset_cts"]]

# parameters
samples <- snakemake@params[['samples']]
sample_id <- snakemake@params[['sample_id']]
ens2geneID <- snakemake@params[['ens2geneID']]

md <- read.delim(samples)

#inputs
counts <- snakemake@input[['counts']]
coverage <- snakemake@input[['coverage']]
st_stats <- snakemake@input[['st_stats']]
meta_data <- md

PCX_counts <- read.delim(counts, row.names=1)
cts <- read.delim(counts, header=TRUE, row.names=1)
cvge <- read.delim(file = coverage)
st_sts <- read.delim(file = st_stats)

#cleaning and pre-processing
counts_ens <- rownames(cts)
counts_ens <- gsub("\\..*","",counts_ens)

cts$gene <- counts_ens
col1 <- colnames(cts)
BR <- col1[1:18]
BR <- gsub('.{21}$',"",BR)
col1 <- gsub("\\_.*","",col1)
colnames(cts) <- col1
colnames(cts)[1:18] <- BR

col2 <- colnames(st_sts)
col2 <- gsub("_bam","",col2)
colnames(st_sts) <- col2
print(colnames(st_sts))

#filtering for selected samples
#pull_sets <- c("CEDAR_REP0158", "BCC_EUS_2019")
#md <- md[md$Set %in% pull_sets,]
cvge <- cvge[cvge$Sample %in% md$SeqID,]
rownames(cvge) <- cvge$Sample
rownames(md) <- md$SeqID
#cts <- data.frame(gene = cts$gene, cts[,colnames(cts) %in% md$SeqID])

md <- md[which(colnames(cts) %in% md$SeqID),]

#putting stats in metadata
#colnames(st_sts)[2:ncol(st_sts)] <- stringr::str_replace(colnames(st_sts)[2:ncol(st_sts)],"_bam","")
#fullst <- st_sts[,1]
#fullst <- as.data.frame(fullst)
#ordered_st <- st_sts[,md$SeqID]
#st_sts <- cbind(fullst, ordered_st)
st_sts[,2:ncol(st_sts)] <- st_sts[,md$SeqID]
#print(colnames(st_sts)[2:ncol(st_sts)] %in% md$SeqID)
print(md$SeqID)
stopifnot(colnames(st_sts)[2:ncol(st_sts)] %in% md$SeqID)
md$input_reads <- t(st_sts[5,2:ncol(st_sts)])
md$input_reads <- as.numeric(md$input_reads)
rownames(md) <- md$SeqID

## Format data ##
cts <- cts[!is.na(cts$gene),]
cts$gene <- as.character(cts$gene)

# Consolidate duplicate gene names
datam <- melt(cts,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
rownames(datac) <- datac[,1]
cts <- datac[,-1]

## Order data the same
cts <- cts[,order(colnames(cts))]
md <- md[order(md$SeqID),]
cvge <- cvge[order(cvge$Sample),]


#check the names are matched up
stopifnot(colnames(cts)==rownames(md) & rownames(md)==rownames(cvge))

# New graph theme
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #panel.border = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(color='black', size=18),
    #axis.text.x=element_text(color='black',size=8),
    panel.grid.major.y = element_line(colour = "black"),
    panel.grid.minor.y = element_line(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_text(angle=45, hjust=1, colour="black"),
    axis.text.y=element_text(colour="black"),
    plot.title=element_text(size=18, face="bold", hjust = 0.5)
  )

dodge <- position_dodge(width = 0.6)
theme_update(plot.title = element_text(hjust = 0.5))

# Import intron and exon fraction results into md table
iv <- match(md$SeqID, cvge$Sample)
md$exon_fraction <- round(cvge[iv,]$Exon, 3)
md$intron_fraction <- round(cvge[iv,]$Intron, 3)
md$intergenic_fraction <- round(cvge[iv,]$Intergenic, 3)

# Now visualize in a different way, retaining the intergenic fraction as well, representing as a stacked barplot
df1 <- dplyr::select(md, intron_fraction)
df1$attribute <- "intronFraction"
names(df1) <- c("Value","attribute")
df1$Sample <- rownames(df1)
df2 <- dplyr::select(md, exon_fraction)
df2$attribute <- "exonFraction"
names(df2) <- c("Value","attribute")
df2$Sample <- rownames(df2)
df3 <- dplyr::select(md, intergenic_fraction)
df3$attribute <- "intergenicFraction"
names(df3) <- c("Value","attribute")
df3$Sample <- rownames(df3)

# Bind these three dataframes together, and then add in other information about each sample for plotting
forPlot <- rbind(df1,df2,df3)
iv <- match(forPlot$Sample, md$SeqID)
head(md[iv,])
forPlot$Status <- md[iv,]$Group
forPlot$CEDAR_ID <- md[iv,]$SeqID
forPlot$Set <- md[iv,]$Set

forPlot$attribute <- factor(forPlot$attribute, levels=c("intergenicFraction","intronFraction","exonFraction"))

fractionBar <- ggplot(forPlot, aes(x=Sample, y=Value, fill=attribute)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each gene attribute") +
  xlab("Sample") +
  scale_fill_brewer( palette = "BuPu", labels = c("Intergenic fraction","Intron fraction","Exon fraction") ) +
  theme(plot.margin = unit(c(0.3,0,0.3,0),"cm"),
        axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

pdf(exon_fraction, 9, 5)
print(fractionBar)
dev.off()

# Barplots of different QC metrics
md$SeqID <- factor(md$SeqID, levels = md$SeqID)
inputReadsBar <- ggplot(data=md, mapping=aes(x=SeqID,y=input_reads,fill=Set)) +
  geom_bar(stat="identity") +
  ylab("Number of input reads") +
  xlab("Sample") +
  theme(plot.margin = unit(c(0.3,0,0.3,0),"cm"),
        axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
inputReadsBar

pdf(input_reads, 9, 5)
print(inputReadsBar)
dev.off()



## Read in raw counts table and make barplot for unique reads lib size
raw <- PCX_counts
raw <- raw[,md$SeqID]
stopifnot(colnames(raw)==md$SeqID)
raw <- raw[,order(colnames(raw))]
md$libSize <- colSums(raw)

uniqueReadsBar <- ggplot(data=md, mapping=aes(x=SeqID,y=libSize,fill=Set)) +
  geom_bar(stat="identity") +
  ylab("Number of gene counts") +
  xlab("Sample") +
  geom_hline(yintercept = 500000, linetype = "dashed") +
  theme(plot.margin = unit(c(0.3,0,0.3,0),"cm"),
        axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
uniqueReadsBar

pdf(unique_reads, 9, 5)
print(uniqueReadsBar)
dev.off()


#library size
md$PASS <- ifelse(md$exon_fraction > 0.34 & md$libSize > 200000, TRUE, FALSE)
pdf(lib_size)
ggplot(data = md, mapping = aes(x = libSize, y = exon_fraction, colour = PASS)) +
  geom_point() +
  geom_vline(xintercept = 200000, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0.34, linetype = "dashed", colour = "grey70")
dev.off()

write.table(md, md_stats, sep = "\t", row.names = F, quote = F)

## PCA plot
genecount <- log(cts+1,2)
# Which values in genecount are greater than 1.5?
thresh <- genecount  > 1.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 5
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- genecount[keep,]
summary(keep)
dim(counts.keep)
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(genecount, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- genecount[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
Genes_wo_hit <- apply(highly_variable_lcpm,1,function(x) all(x==0))
genecount_nozero <- highly_variable_lcpm[!Genes_wo_hit,]
genecount_PCA <- t(genecount_nozero)

# Run PCA
PCA_out <- prcomp(genecount_PCA, scale. = TRUE)

Groups <- md$Group
# Generate the summary stats for prcomp
eigs <- PCA_out$sdev^2
summary <- rbind(
  SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))
# Find the proprtion of PC1 and PC2
propPC1 <- round(as.numeric(summary[2,1]) * 100)
propPC2 <- round(as.numeric(summary[2,2]) * 100)

# Plot PCA
my.pca <- data.frame(PCA_out$x[,1:3])
my.pca$Sample <- rownames(my.pca)
iv <- match(my.pca$Sample, md$SeqID)
potential_batch <- c("Set","input_reads","libSize","exon_fraction","PASS","Group")
my.pca <- data.frame(my.pca, md[iv,potential_batch])

ggplot(data = my.pca, mapping = aes(x = PC1,y = PC2, colour = Group)) +
  geom_point(size = 2) +
  theme_bw() +
  xlab(paste0("PC-1 (", propPC1, "%)")) +
  ylab(paste0("PC-2 (", propPC2, "%)"))

## Now let's label by library size and by PASS

plots <- list()
for (i in 1:length(potential_batch)) {
  if (i==1) {
    p <- ggplot(data = my.pca, mapping = aes(x = PC1,y = PC2, colour = Set)) +
      geom_point(size = 3) +
      theme_bw() +
      xlab(paste0("PC-1 (", propPC1, "%)")) +
      ylab(paste0("PC-2 (", propPC2, "%)"))
  } else {
    p <- ggplot(data = my.pca, mapping = aes_string(x = "PC1",y = "PC2", colour = potential_batch[i])) +
      geom_point(size = 3)  +
      theme_bw() +
      xlab(paste0("PC-1 (", propPC1, "%)")) +
      ylab(paste0("PC-2 (", propPC2, "%)"))
  }
  plots[[potential_batch[i]]] <- p
}

pdf(batch_pca, width = 22, height = 12)
plot_grid(plotlist = plots, labels = potential_batch)
dev.off()

tmp <- sweep(cts, 2, colSums(raw), '/') * 1e6
PCX_RPM <- tmp
write.csv(PCX_RPM, sub_counts)
sub_cts <- cts[,md$SeqID[which(md$PASS == TRUE)]]
sub_cts <- round(sub_cts)
write.table(sub_cts, subset_cts)
