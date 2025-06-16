rule sort_by_coord:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam"
    output:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sortByCoord.rmd.bam"
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        "samtools sort -O bam {input} -o {output}"

rule index_rmd_sortByCoord:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sortByCoord.rmd.bam"
    output:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sortByCoord.rmd.bam.bai"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "samtools index {input} {output}"

rule bigWig:
    input:
        bam = "samples/genecounts_rmdp/{sample}_bam/{sample}_sortByCoord.rmd.bam",
        bai = "samples/genecounts_rmdp/{sample}_bam/{sample}_sortByCoord.rmd.bam.bai"
    output:
        "samples/bigwigs/{sample}.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output}"
