contrast=get_contrast
rule debks_prep_gtf:
    input:
        gtf = config["gtf_file"]
    output:
        "tools/sort.gtf.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        """grep -v '#' {input} |sort -k 1,1 -k 4,4n |bgzip > tools/sort.gtf.gz
        tabix sort.gtf.gz"""

rule debks_prep_genome:
    input:
        genome=config["genome"]
    output:
        "{genome}.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """samtools faidx {input}"""

rule contrast_namer:
    input:
        expand("ciri2_output/{contrast}/ciri2out_bwa_map_{sample}.txt", sample=SAMPLES, contrast=config["diffexp"]["contrasts"])
    output:
        "ciri2_output/{contrast}/{contrast}_namelist.txt",
        "ciri2_output/{contrast}/{contrast}_namelengths.txt"
    params:
        samples=config["omic_meta_data"],
        sample_id=config["sample_id"],
        linear_model=config["linear_model"],
        contrast=get_contrast
    script:
        "../scripts/name_lister.R"

rule debks_merge:
    input:
        "ciri2_output/{contrast}/{contrast}_namelist.txt"
    output:
        "results/debks_output/{contrast}/merge_pos.txt",
        "results/debks_output/{contrast}/merge_linear.txt",
        "results/debks_output/{contrast}/merge_circ.txt"
    params:
        contrast=config["diffexp"]["contrasts"]
    conda:
        "../envs/debks.yaml"
    shell:
        """DEBKS merge -s CIRI2 -f {input} -o results/debks_output/{wildcards.contrast}/merge"""

rule debks_anno:
    input:
        pos = "results/debks_output/{contrast}/merge_pos.txt",
        genome = config["genome"],
        sort_gtf = "tools/sort.gtf.gz"
    output:
        "results/debks_output/{contrast}/anno_len.txt"
    params:
        contrast=get_contrast
    conda:
        "../envs/debks.yaml"
    shell:
        """DEBKS anno -c {input.pos} -m {input.genome} -g {input.sort_gtf} -o results/debks_output/{wildcards.contrast}/anno"""

rule debks_dec:
    input:
        circ = "results/debks_output/{contrast}/merge_circ.txt",
        lin = "results/debks_output/{contrast}/merge_linear.txt",
        lengths = "ciri2_output/{contrast}/{contrast}_namelengths.txt",
        anno = "results/debks_output/{contrast}/anno_len.txt"
    output:
        "results/debks_output/{contrast}/{contrast}_dec_circRNA.txt"
    params:
        contrast=get_contrast
    conda:
        "../envs/debks.yaml"
    shell:
        """lens=$(head -n 1 {input.lengths})
        DEBKS dec -c {input.circ} -l {input.lin} -n $lens -f 1 -e {input.anno} -o {output}"""
