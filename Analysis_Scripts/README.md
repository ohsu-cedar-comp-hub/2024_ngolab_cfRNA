# Reproduction of paper methods and figures

Use **main_script.R** to run through all data analysis and figure generation. This will step through the normalization of the raw count matrices, DE gene discovery, classification tasks, tissue deconvolution of DE gene signals, surival analysis, and generation of all main and supplemental figures presented in the paper.

We also provide **no_intrinsic_norm.R**, which runs through a similar biomarker discover and classification analysis using a more traditional TPM normalization.

These scripts should be run in the Analysis_Scripts folder, as they will depend on the R files and data files there. A breakdown of the data files included is provided below.

# Cf-normalization
Our cf-normalization method can be used outside of this analysis by running the cf_norm() function in normalization.R. This function takes a raw count matrix of reads as the first argument (genes in rows, samples in columns), and a data frame of intrinsic and extrinsic factors as the second argument. We provide these factors as learned on a healthy cohort in the file spin_data_intrinsic.csv. Alternatively the factors in updated_intrinsic.csv can also be used. This file has additional factor estimates for genes that were present in the PDAC data but not in the healthy reference data.

See main_script.R for an example use of cf_norm and loading the factor files.

# Data files provided
pdac_genecount.csv.gz: Zipped raw gene count matrix for all CEDAR and BCC samples analyzed.

pdac_meta.csv: Metadata file for pdac_genecount.csv.gz.

PDAC_survival.csv: Additional metadata indicating the disease status and survival of PDAC patients at followup date.

pdac_qc_metrics.csv: QC metrics calculated for RNA processing of CEDAR and BCC samples.

normalized_pdac_genecount.csv: Normalized gene counts (using cf-normalization) for all CEDAR and BCC samples. Uses same metadata file as pdac_genecount.csv.gz.

reference_data_counts.csv: Reference dataset of healthy donors, with plasma undergoing 10 different processing and handling conditions.  Used to estimate intrinsic and extrinsic gene factors for normalization.

reference_metadata.csv: Metadata for reference dataset. Includes donor ID code and handling conditions for each sample.

spin_data_instrinsic.csv: Estimated intrinsic and extrinsic factors computed for each gene in reference cohort.

updated_intrinsic.csv: Estimated intrinsic and extrinsic factors computed for each gene in reference cohort, plus additional estimates for genes only present in pdac cohort. 

Tissue_RNA_Atlas.csv.gz: Zipped matrix of TPM normalized tissue RNA values. Used for DE gene source and deconvolution analysis. Downloaded from “RNA HPA tissue gene data” from Human Protein Atlas version 22.0

hsapiens_gene_ensemble__gene__main.txt.gz: Zipped gene ensemble file (hg38).

g_to_e.csv: Table of gene name, ensemble name, and transcipt length. Used for TPM normalization.
