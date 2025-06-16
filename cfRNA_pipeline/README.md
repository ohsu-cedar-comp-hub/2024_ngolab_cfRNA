# cf-RNA Processing Pipeline

## Setup

To generate RNA counts from the raw fastq sequencing files we used this Snakemake pipeline that was originally worked on by Breeshey Roskams-Hieter. Several pre-requisite files and pieces of software are necessary for running the pipeline, and should be included in the omic_config.yaml file:

- A GTF file (gencode version 27 used here)
- A genome file (UCSC hg38 used here)
- A bed file of exon locations for QC
- A star index of the genome used
- A GTF of housekeeping gene exons
- fastqscreen
- STAR aligner (ver. 2.5.3a used here)
- sickle tool
- picard tool, pointing to the MarkDuplicates.jar
- bamstats tool
 
Once these locations are added to the omic_config.yaml file, conda should be set up if it hasn't already been installed on the system. The pipeline was set up to run from a snakemake conda environment, and submit the processing jobs in a slurm workflow manager.
 
After this software setup is complete, the final step before running the pipeline is to link the raw, gzipped fastq files into a sub-directory located in ./samples/raw/
Each soft-linked sample should have a sample name, followed by a "_" and either a "R1" or an "R2" e.g. a paired end sample named 001 would have two soft-links in the ./samples/raw/ directory named "001_R1.fq" and "001_R2.fq" that would point to the gzipped R1 and R2 fastq files for that sample.

### Linking example

Symlinking allows you to run the pipeline without having to move files around and copy them.
It can be performed as demonstrated below for each file.

```{bash}

$ ln -s /source/to/main_file_1.fq.gz /new/file/here_R1.fq

```
