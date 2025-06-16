library("dplyr")
library("DESeq2")
library("data.table")

counts = snakemake@input[['counts']]

metadata <- snakemake@params[['samples']]

sampleID <- snakemake@params[['sample_id']]

Type <- snakemake@params[['linear_model']]

contrast <- snakemake@params[['contrast']]

baseline <- contrast[[2]]

target <- contrast[[1]]

output = snakemake@output[['rds']]

rld_out = snakemake@output[['rld_out']]

list_output = snakemake@output[['list_rds']]

list_rld_out = snakemake@output[['list_rld_out']]

ERCC_output = snakemake@output[['ERCC_rds']]

ERCC_rld_out = snakemake@output[['ERCC_rld_out']]

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Read in metadata table and order according to sampleID
md <- read.delim(file=metadata, stringsAsFactors = FALSE)
md <- md[order(md[sampleID]),]

# Read in counts table
subdata <- read.table(counts, row.names=1)
counts_ens <- rownames(subdata)
rownames(subdata) <- gsub("\\..*","",counts_ens)
subdata <- subdata[,order(colnames(subdata))]

# Extract only the Types that we want in further analysis & only the PP_ID and Status informative columns
md <- filter(md, !!as.name(Type) == baseline | !!as.name(Type) == target, !!as.name(sampleID) %in% colnames(subdata))

# Keep only the PP_IDs of the types we have chosen in the metadata table above
rownames(md) <- md[[sampleID]]
md[[sampleID]] <- NULL
keep <- colnames(subdata)[colnames(subdata) %in% rownames(md)]
subdata <- subdata[, keep]
dim(subdata)
subdata <- as.matrix(subdata)
# Check
stopifnot(rownames(md)==colnames(subdata))

# Obtain the number of genes that meet padj<0.01 for reference line in histogram
dds <- DESeqDataSetFromMatrix(countData=subdata,
                              colData=md,
                              design= as.formula(paste('~',Type)))

dds <- estimateSizeFactors(dds)

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Normalization and pre-processing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=output, version = 2)

# obtain normalized counts
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=rld_out, version = 2)

## Normalizing to gene List
dds <- DESeqDataSetFromMatrix(countData=subdata,
                              colData=md,
                              design= as.formula(paste('~',Type)))
genelist <- c("ENSG00000198888", "ENSG00000198763", "ENSG00000198899", "ENSG00000274012", "ENSG00000161016", "ENSG00000244734", "ENSG00000034510", "ENSG00000108298", "ENSG00000142541", "ENSG00000167526", "ENSG00000142937", "ENSG00000177954", "ENSG00000168298", "ENSG00000177600", "ENSG00000205542", "ENSG00000115268", "ENSG00000188536", "ENSG00000149806", "ENSG00000133112", "ENSG00000140988")
genenums <- which(rownames(counts(dds)) %in% genelist)
dds <- estimateSizeFactors(dds, controlGenes = genenums)

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Normalization and pre-processing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=list_output, version=2)

# obtain normalized counts
list_rld <- rlog(dds, blind=FALSE)
saveRDS(list_rld, file=list_rld_out, version=2)

## Normalizing to ERCC counts
dds <- DESeqDataSetFromMatrix(countData=subdata,
                              colData=md,
                              design= as.formula(paste('~',Type)))

ERCC_rows <- which(rownames(dds) %like% "ERCC")
dds <- estimateSizeFactors(dds, controlGenes = ERCC_rows)

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Normalization and pre-processing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=ERCC_output, version=2)

# obtain normalized counts
ERCC_rld <- rlog(dds, blind=FALSE)
saveRDS(ERCC_rld, file=ERCC_rld_out, version=2)

