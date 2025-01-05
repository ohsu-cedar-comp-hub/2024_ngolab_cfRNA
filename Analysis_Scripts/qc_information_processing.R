library(tidyverse)

#changed one
#change these filepaths to be the location of where the following files are on your computer
####Loading QC metadata
qc_metrics_metadata <- read.csv("~/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/PDAC_paper_scripts/PDAC_paper_meta.csv")

#read.csv(./qcmetadata)
#Filter by the samples that we are using
qc_metrics <- read.delim("~/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/PDAC_analysis/ncRNA_PDAC_Metadata_W_QC_Metrics.txt") %>%
  as.data.frame()
qc_metrics_filtered <- qc_metrics %>% filter(SeqID %in% qc_metrics_metadata$SeqID)
qc_metrics_filtered$Set <-ifelse(grepl("CEDAR", qc_metrics_filtered$Set), "CEDAR", "BCC")
qc_metrics_filtered <- qc_metrics_filtered # #only if we want to make it waterfall%>% arrange(Set, -libSize)
qc_metrics_filtered_correct_col <- qc_metrics_filtered %>% 
  select(-c("SampleID", "Clinical.Diagnosis")) 
qc_metrics_filtered$SeqID <- factor(qc_metrics_filtered$SeqID)
#write.table(qc_metrics_filtered, file=".")

#Loading Unprocessed sequence counts
#need to get the counts that aren't just the passing counts
#pdac_data_unprocessed<- read.csv("~/Downloads/paper_norm_counts.csv", sep=",", row.names = 1)
pdac_counts <- read.csv("~/Downloads/PDAC_5set_genecount.csv")
pdac_counts <- pdac_counts %>% pivot_longer(!gene_name) %>%
  group_by(gene_name,name) %>% summarize(value=sum(value)) %>% pivot_wider()
pdac_data_unprocessed <- pdac_counts %>% as.data.frame() %>% filter(!is.na(gene_name))
rownames(pdac_data_unprocessed) <- pdac_data_unprocessed$gene_name
pdac_data_unprocessed$gene_name <- NULL
#gene information about ensembl 90 from stored ensembl database
gene_info_ens90<- read.delim("./hsapiens_gene_ensembl__gene__main.txt.gz", header = F)

#########################
## Creating fraction table
#TOKEEP
qc_fraction_info <- qc_metrics_filtered %>% 
  select(c("exon_fraction","intron_fraction", "intergenic_fraction", "libSize",
           "input_reads", "SeqID")) %>% pivot_longer(!SeqID)
fraction_info_summary <- qc_fraction_info %>% group_by(name) %>%
  summarize(mean=mean(value), max=max(value), min=min(value))
fraction_info_summary$name <- c("Exon Fraction", "Input Reads", "Intergenic Fraction","Intron Fraction", "Library Size")
#write.table(fraction_info_summary, file="./PDAC_read_type_fraction.txt", )


#creating plot of total reads and deduplicated
qc_size<- qc_metrics_filtered %>%
   dplyr::select(c(libSize, SeqID, Set)) 
qc_size_longer<- pivot_longer(qc_size, !c("SeqID", "Set"))
qc_size_longer$Set <- factor(qc_size_longer$Set)
qc_size_longer$SeqID <- factor(qc_size_longer$SeqID, levels =qc_metrics_filtered$SeqID )
color_vec <- c('#66c2a5','#8da0cb')
names(color_vec) <- c("CEDAR", "BCC")
#now we are plotting the qc curves for both of the cohorts
qc_size_plot_CEDAR<- ggplot(qc_size_longer %>% filter(Set=="CEDAR"), aes(x=SeqID, y=value, fill=Set)) + 
  geom_bar(stat="identity") + theme_classic() +xlab("Sample")+
  ylab("Total deduplicated reads") +
  theme(axis.text.x = element_text(angle = 90, size=10), legend.position = "none")+
  ggtitle("CEDAR Samples")+
  scale_fill_manual(values=color_vec)+
  ylim(0,8000000)
qc_size_plot_BCC <- ggplot(qc_size_longer%>%filter(Set=="BCC"), aes(x=SeqID, y=value, fill=Set)) + 
  geom_bar(stat="identity") + theme_classic() +xlab("Sample")+
  ylab("Total deduplicated reads") +
  theme(axis.text.x = element_text(angle = 90, size=10), legend.position = "none")+
  ggtitle("BCC Samples")+
  scale_fill_manual(values=color_vec)+
  ylim(0,8000000)
#want to make it easy to see any differences without being confused by the scale^

## Creating intron exon intergenic fraction plot
#this plot will be used with the saturation curves to show
qc_intron_exon_fraction <- qc_metrics_filtered %>%
   dplyr::select(c(intron_fraction, exon_fraction, intergenic_fraction, SeqID, Set)) 
qc_fractions <- pivot_longer(qc_intron_exon_fraction, !c("SeqID", "Set"))                    
qc_fractions$SeqID <- factor(qc_fractions$SeqID, levels=qc_metrics_filtered$SeqID)
qc_fractions$name <- factor(qc_fractions$name, levels=c("intergenic_fraction","intron_fraction","exon_fraction"))
color_list <- c('#66c2a5','#fc8d62','#8da0cb')
names(color_list) <- c("intergenic_fraction", "exon_fraction", "intron_fraction")
qc_fractions_plot_BCC <- ggplot(qc_fractions %>% filter(Set=="BCC"), aes(x=SeqID, y=value, fill=name)) + 
  geom_bar(stat="identity", position="fill") +theme_classic() + scale_fill_manual(values=color_list)+
  theme(axis.text.x = element_text(angle = 90, size=10)) +
  ylab("Fraction of total reads") +xlab("")
qc_fractions_plot_CEDAR <- ggplot(qc_fractions%>%filter(Set=="CEDAR"), aes(x=SeqID, y=value, fill=name)) + 
  geom_bar(stat="identity", position="fill") +theme_classic() + scale_fill_manual(values=color_list)+
  theme(axis.text.x = element_text(angle = 90, size=10)) +
  ylab("Fraction of total reads") +xlab("")

##################################
##################################
##Working with mrna Sequences
library(vegan)
#Take out and remove the ERCC transcripts that are present, we will not be using 
cohorts <- qc_metrics_metadata$SeqID #select only the samples included in the paper analysis
counts <- t(pdac_data_unprocessed[, cohorts])
#Remove columns of genes where all of the counts are 0
counts <- counts[, colSums(counts) > 0]

###########################################
##Protein Coding and biotype by sample graph
counts_df <- counts %>% as.data.frame()
genes <- colnames(counts_df)
counts_df$sample_name <- rownames(counts_df)
long_counts <- pivot_longer(counts_df, !sample_name)

#we use gencode v 27 which corresponds with ensembl 90
#get version 90 from downloaded ensembl information for genes main from database online
#downloaded from https://ftp.ensembl.org/pub/release-90/mysql/ensembl_mart_90/
biotype <- gene_info_ens90[match(long_counts$name, gene_info_ens90$V8), "V3"]
long_counts$biotype <- biotype

#here we check to see what types are the most common
total_amounts<- long_counts %>% group_by(biotype) %>%
  summarize(total_counts=sum(value))
#it looks like most reads fall within one of the following four categories
long_counts$biotype <- ifelse(long_counts$biotype %in% c("protein_coding",
                                                      "Mt_rRNA",
                                                      "processed_pseudogene",
                                                      "lincRNA", 
                                                      "misc_RNA"),
                                                      long_counts$biotype, "other")
#reorder the factor levels so that they are in order of most to least prevalent
ordered_biotypes = rev(c("protein_coding",
"Mt_rRNA",
"processed_pseudogene",
"lincRNA", 
"misc_RNA", "other"))
categorical_palette <- c('#8da0cb','#a6d854', '#e78ac3','#66c2a5','#ffd92f','#fc8d62')
long_counts$biotype <- factor(long_counts$biotype, levels= ordered_biotypes)
long_counts_grouped <- long_counts %>% group_by(biotype, sample_name) %>%
  summarize(total_counts=sum(value))
#now plot a bar plot for protein coding
long_counts_grouped$group <- qc_metrics_filtered[match(long_counts_grouped$sample_name, qc_metrics_filtered$SeqID), "Group"]
long_counts_grouped$set <- qc_metrics_filtered[match(long_counts_grouped$sample_name, qc_metrics_filtered$SeqID), "Set"]
long_counts_grouped$sample_name <- factor(long_counts_grouped$sample_name, levels=qc_metrics_filtered$SeqID)

#creating separate plots for each cohort
biotype_qc_plot_BCC <- ggplot(long_counts_grouped %>%filter(set=="BCC") , aes(x=sample_name, y=total_counts, fill=biotype)) +
  geom_bar(stat="identity", position="fill") +
  theme_classic()+
  scale_fill_manual(values=categorical_palette)+
  theme(axis.text.x = element_text(angle = 90, size=10)) + 
  ylab("Fraction of total reads") +xlab("")+
  theme_classic()
biotype_qc_plot_CEDAR <- ggplot(long_counts_grouped%>%filter(set=="CEDAR"), aes(x=sample_name, y=total_counts, fill=biotype)) +
  geom_bar(stat="identity", position="fill") +
  theme_classic()+
  scale_fill_manual(values=categorical_palette)+
  theme(axis.text.x = element_text(angle = 90, size=10)) + 
  ylab("Fraction of total reads") +xlab("")
#lets create our individual by sample qc plots that combine everything
qc_size_bar_CEDAR <- qc_size_plot_CEDAR +coord_flip() #+theme(legend.position="top")
qc_size_bar_BCC <- qc_size_plot_BCC +coord_flip() #+theme(legend.position="top")
qc_fraction_bar_CEDAR <- qc_fractions_plot_CEDAR + coord_flip()
qc_fraction_bar_BCC <- qc_fractions_plot_BCC + coord_flip()
qc_biotype_bar_CEDAR <- biotype_qc_plot_CEDAR + coord_flip()
qc_biotype_bar_BCC <- biotype_qc_plot_BCC+ coord_flip()
library(cowplot)

plot_grid(
  qc_size_bar_CEDAR,
  qc_fraction_bar_CEDAR,
  qc_biotype_bar_CEDAR,
  labels = c('A', 'B', 'C'),
  align="h",
  nrow=1,
  ncol=3
) 
plot_grid(
  qc_size_bar_BCC,
  qc_fraction_bar_BCC,
  qc_biotype_bar_BCC,
  labels = c('A', 'B', 'C'),
  align="h",
  nrow=1,
  ncol=3
)

#now compare biotypes for broader groups and done with just per sample plotting
#comparing percentages across sample and sample types 
percent_long_counts <-long_counts 
percent_long_counts<- percent_long_counts %>% group_by(sample_name) %>%
  mutate(value=value/sum(value)) %>% group_by(biotype, sample_name) %>% 
  summarize(total_counts= sum(value))
#get the diagnosis/group information
percent_long_counts$group <- qc_metrics_metadata[match(percent_long_counts$sample_name, qc_metrics_filtered$SeqID), "Group"]
#get the source and rename it
percent_long_counts$cohort <- qc_metrics_metadata[match(percent_long_counts$sample_name, qc_metrics_filtered$SeqID), "Source"]
percent_long_counts$cohort <- ifelse(percent_long_counts$cohort=="CEDAR_2020", "CEDAR", "BCC")
#lets also plot it for the discovery vs validation cohort
#TOKEEP
biotype_by_cohort_group <- ggplot(percent_long_counts, aes(x=group, y=total_counts, fill=biotype, color=biotype)) +
  geom_boxplot() + facet_wrap(~cohort, ncol=1) +
  ylab("% of total read counts in sample") +xlab("Sample Type") + theme_minimal()

#write percent biotype table
summarized_percent_biotype <- percent_long_counts %>% group_by(biotype) %>%
  summarize(minimum=min(total_counts), mean=mean(total_counts), maximum=max(total_counts))
write.table(summarized_percent_biotype, file="./PDAC_biotype_information_table.txt")

summarized_percent_biotype_sample<- long_counts%>% group_by(sample_name) %>%
  mutate(value=value/sum(value)) %>% group_by(biotype, sample_name) %>% 
  summarize(total_counts= sum(value)) %>% 
  pivot_wider(names_from=biotype, values_from=total_counts)
#now lets create a total datatable for this information
combined_table <- merge(summarized_percent_biotype_sample, qc_metrics_filtered,by.y="SeqID", by.x="sample_name" ) %>%
  select(c("other", "misc_RNA", "lincRNA", "processed_pseudogene",
           "Mt_rRNA", "protein_coding", "exon_fraction", "intergenic_fraction",
           "intron_fraction", "libSize", "sample_name", "Set" ))
colnames(combined_table) <- c("other biotype", "misc_RNA", "lincRNA", "processed_pseudogene",
                               "Mt_rRNA", "protein_coding", "Exon fraction", "Intergenic fraction",
                               "Intron fraction", "Total deduplicated reads", "Sample name", "Sample source" )

combined_table <- combined_table[,c("Sample name", "Sample source", "Total deduplicated reads",
"Exon fraction", "Intergenic fraction","Intron fraction",
"protein_coding", "Mt_rRNA",
"processed_pseudogene",  "lincRNA",
"misc_RNA", "other biotype" )]

combined_table <- combined_table %>% arrange(`Sample source`, `Sample name`)
###################################
##Saturation Curves
#this function is a very lightly modified version of the vegan::rarecurve function
#This function from vegan simply creates a rarefaction curve. Some changes were needed
# in order for this function to be usable with multiple different samples in a way that was easy
#To graph
rarecurve_t <- function(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
                        label = TRUE, col, lty, ...){
  ## matrix is faster than data.frame
  x <- as.matrix(x)
  ## check input data: must be counts
  if (!identical(all.equal(x, round(x)), TRUE))
    stop("function accepts only integers (counts)")
  ## sort out col and lty
  if (missing(col))
    col <- par("col")
  if (missing(lty))
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  ## remove empty rows or we fail
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0,, drop =FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  ## rep col and lty to appropriate length
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  ## Rarefy
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      ## don't want names on n an `c` adds a name from `tot[i]`)
      n <- c(n, tot[i], use.names = FALSE)
    }
    drop(rarefy(x[i,], n))
  })
  #lets make sure that our names are still on there
  names(out) <- rownames(x)
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  
  return(out)
  
  
}

#make the counts round to a whole number and remove fractional counts
#processing saturation curves
sample <- rowSums(counts)
curveout<-rarecurve_t(counts, step=10000, sample)

#process the weird format that the names come out of the rarecurve_t function
#Making sure to get the names of the specific samples that we care about
get_dataframe <- function(curveout_result, name){
  names_result <- names(curveout_result) %>% stringr::str_sub(., start=2)
  result<- as.vector(unlist(curveout_result))
  df <- data.frame(rarefaction=names_result, V2=result)
  colnames(df)[names(df)=="V2"] <- name
  return(df)
}

#Look through all of the different samples and pull out this function to create a list of dataframes
dataframes <- lapply(names(curveout),function(x) get_dataframe(curveout[[x]],x))
#combind all of these dataframes together This is slightly annoying because it has to be merge because of the different
#Row length that the different samples have
reduced <- purrr::reduce(dataframes, function(x,y) merge(x,y, by.x="rarefaction", by.y="rarefaction", all.x=T, all.y=T))
#Make sure that everything is all numbers
reduced$rarefaction <- as.numeric(reduced$rarefaction)
reduced <- pivot_longer(reduced, cols= !contains("rarefaction"))
reduced <- reduced[!is.na(reduced$value),]
#lets add another column that will make things fractional
reduced_w_fractional<- reduced %>% group_by(name) %>% 
  mutate(fractional_rarefaction=rarefaction/max(rarefaction)) %>%
  ungroup()
reduced <- reduced_w_fractional
#Plot Saturation Curves
type <- match(reduced$name, qc_metrics_metadata$SeqID)
reduced$type <- qc_metrics_metadata[type, c("Group")] %>% str_replace("Benign pancreas", "Benign Pancreas")
reduced <- reduced[order(reduced$type), ]
reduced$name <- factor(reduced$name, levels= unique(reduced$name))
color_palette_base <- c('#F6C142','#C2D7EC','#4075B1', '#4075B1', '#DF8244','#B02318', '#68379a', '#4eac5b')
names(color_palette_base) <- c("IPMN", "Acute pancreatitis","Chronic pancreatitis", "Pancreatitis",
                               "Cancer_other", "Islet Cell Tumor","PDAC", "Benign Pancreas")
rarefaction_curve_overall <- ggplot(data=reduced, aes(x=fractional_rarefaction, y=value,group=name, color=type)) +
  geom_line(show.legend=T, alpha=0.5) +
   scale_color_manual(values=color_palette_base) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))
#labelling and facet wrapping options
#+facet_wrap(~name)
#geom_label(data= reduced%>% group_by(.,name)%>%filter(value==max(value)) %>% ungroup(), aes(label=name), show.legend = F) +
#TOKEEP
rarefaction_curve_plot <- rarefaction_curve_overall +ylab("Total Number of Unique Genes In Sample") +xlab("Fraction of Reads Sampled")


