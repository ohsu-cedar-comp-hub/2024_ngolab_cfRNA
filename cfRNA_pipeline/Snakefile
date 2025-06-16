__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")
print(SAMPLES)

#removed the index col for our metadata no longer needed it
sample_id = config["sample_id"]
md = pd.read_table(config["omic_meta_data"], index_col = sample_id, dtype=str)
condition = config["linear_model"]
baseline = config["TE_baseline"]

# TO FILTER 
control = md[md[condition]==baseline]
treat = md.loc[~md.index.isin(control.index)] 

control_paths = ['samples/star_TE/{}/Aligned.out.bam'.format(x) for x in control.index]

treat['cond_paths'] = ['samples/star_TE/{}/Aligned.out.bam'.format(x) for x in treat.index]
cond_paths = {key:x['cond_paths'].values.tolist() for key,x in treat.groupby(condition)}
CONDITIONS = cond_paths.keys()

# Wildcard function to grab proper condition
def get_TE(wildcards):
    return cond_paths[wildcards.condition]

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']


with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)


result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(config["omic_meta_data"]) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]

full_contrast = config["diffexp"]["contrasts"]

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

for condition in CONDITIONS:
    message("Condition " + condition + " will be processed")
print(control_paths)
print(cond_paths)
print(SAMPLES)


rule all:
    input:
#        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_{fastq_ext}_t_fastqc.zip", sample = SAMPLES, fastq_ext = fastq_ext),
#        expand("samples/fastqscreen/{sample}/{sample}_{fastq_ext}_t_screen.{fastqscreen_ext}", sample=SAMPLES, fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
#        "data/{project_id}_counts_w_stats.txt".format(project_id=config['project_id']),
#        "data/{project_id}_counts.txt".format(project_id=config['project_id']),
#        "data/{project_id}_exon_counts.txt".format(project_id = config["project_id"]),
#        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
#        expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
#        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
#        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
#        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
#        "results/tables/read_coverage.txt",
#        expand("{project_id}_QC/Exon.Fraction.pdf", project_id=config["project_id"]),
#        expand("{project_id}_QC/Input.Reads.Barplot.pdf", project_id=config["project_id"]),
#        expand("{project_id}_QC/Unique.Reads.Barplot.pdf", project_id=config["project_id"]),
#        expand("{project_id}_QC/Scatter.Plot.LibSize.ExonFraction.pdf", project_id=config["project_id"]),
#        expand("{project_id}_QC/{project_id}_Metadata_W_QC_Metrics.txt", project_id=config["project_id"]),
#        expand("{project_id}_QC/{project_id}_PCA_potential_batch_grid.pdf", project_id=config["project_id"]),
#        expand("{project_id}_QC/{project_id}_subset_counts_RPM.txt", project_id=config["project_id"]),
#        expand("{project_id}_QC/{project_id}_subset_counts_pass.txt", project_id=config["project_id"]),
#        expand("samples/star_TE/{sample}/Aligned.out.bam", sample = SAMPLES),
#        expand("results/TEtranscripts/{condition}.cntTable", condition = CONDITIONS),
#        expand("results/TEtranscripts/{condition}_sigdiff_gene_TE.txt", condition = CONDITIONS)
#        expand("ciri2_output/ciri2out_bwa_map_{sample}.txt", sample = SAMPLES),
#        expand("samples/circ_explorer/{sample}/{sample}_known_circ.txt", sample = SAMPLES)
#        "tools/sort.gtf.gz",
#        expand("{genome}.fai", genome = config["genome"]),
#        expand("ciri2_output/{contrast}/{contrast}_namelist.txt", contrast = config["diffexp"]["contrasts"]),
#        expand("ciri2_output/{contrast}/{contrast}_namelengths.txt", contrast = config["diffexp"]["contrasts"]),
#        expand("results/debks_output/{contrast}/{contrast}_merge_pos.txt", contrast=config["diffexp"]["contrasts"]),
#        expand("results/debks_output/{contrast}/{contrast}_merge_linear.txt", contrast=config["diffexp"]["contrasts"]),
#        expand("results/debks_output/{contrast}/{contrast}_merge_circ.txt", contrast = config["diffexp"]["contrasts"]),
#        expand("results/debks_output/{contrast}/{contrast}_dec_circRNA.txt", contrast = config["diffexp"]["contrasts"])
#        "results/debks_output/{contrast}/{contrast}_dec_circRNA.txt".format(contrast=config["diffexp"]["contrasts"])
#        "results/diffexp/group/LRT_pca.pdf",
#        "results/diffexp/group/LRT_pca_listnorm.pdf",
#        "results/diffexp/group/LRT_pca_ERCCnorm.pdf",
#        "results/diffexp/group/MDS_table.txt",
#        "results/diffexp/group/MDS_table_listnorm.txt",
#        "results/diffexp/group/MDS_table_ERCCnorm.txt",
#        "results/diffexp/group/LRT_density_plot.pdf",
#        "results/diffexp/group/LRT_density_plot_listnorm.pdf",
#        "results/diffexp/group/LRT_density_plot_ERCCnorm.pdf",
#        expand(["results/diffexp/pairwise/{contrast}.qplot.pdf","results/diffexp/pairwise/{contrast}.qhist.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv"],contrast=config["diffexp"]["contrasts"]),
#        expand(["results/diffexp/pairwise/{contrast}.qplot.listnorm.pdf","results/diffexp/pairwise/{contrast}.qhist.listnorm.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.listnorm.tsv"],contrast=config["diffexp"]["contrasts"]),
#        expand(["results/diffexp/pairwise/{contrast}.qplot.ERCCnorm.pdf","results/diffexp/pairwise/{contrast}.qhist.ERCCnorm.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.ERCCnorm.tsv"],contrast=config["diffexp"]["contrasts"]),
#        expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
#        expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO_listnorm.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO_listnorm.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
#        expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO_ERCCnorm.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO_ERCCnorm.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
#        expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
#        expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.listnorm.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
#        expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.ERCCnorm.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
#        expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf", contrast = config["diffexp"]["contrasts"]),
#        expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.listnorm.pdf", contrast = config["diffexp"]["contrasts"]),
#        expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.ERCCnorm.pdf", contrast = config["diffexp"]["contrasts"]),
#        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot.html"],contrast = config["diffexp"]["contrasts"]),
#        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot_listnorm.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot_listnorm.html"],contrast = config["diffexp"]["contrasts"]),
#        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot_ERCCnorm.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot_ERCCnorm.html"],contrast = config["diffexp"]["contrasts"]),
#        "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
#        "results/diffexp/glimma-plots/{project_id}.mds_plot_listnorm.html".format(project_id=project_id),
#        "results/diffexp/glimma-plots/{project_id}.mds_plot_ERCCnorm.html".format(project_id=project_id)
include: "rules/align_rmdp.smk"
#include: "rules/omic_qc.smk"
#include: "rules/QC_init.smk"
#include: "rules/TE.smk"
#include: "rules/circ.smk"
#include: "rules/debks.smk"
#include: "rules/deseq.smk"
