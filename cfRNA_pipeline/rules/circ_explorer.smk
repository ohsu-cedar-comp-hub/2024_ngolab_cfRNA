contrast=get_contrast
rule circ_explorer:
        fwd = "samples/raw/{sample}_R1.fq",
        rev = "samples/raw/{sample}_R2.fq"
    output:
        "circ_explorer_test"
    conda:
        "../envs/circexplorer.yaml"
    shell:
        """fast_circ.py annotate -r /home/groups/CEDAR/spilioto/tools/hg38_ref_all.txt -g /home/groups/CEDAR/spilioto/Projects/PDAC_ncRNA/tools/GRCh38.primary_assembly.genome.fa -G /home/groups/CEDAR/spilioto/tools/hg38_ref_all.gtf --pe {rev} -o circ_explorer_test -f {fwd}"""
