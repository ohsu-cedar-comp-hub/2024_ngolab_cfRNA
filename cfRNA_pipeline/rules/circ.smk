rule bwa: 
        input: "samples/raw/{sample}_R1.fq" 
        output:"data/sam/bwa_map_{sample}.sam"  

        params: name="bwa_{sample}", partition="exacloud", mem="64000" 
        threads: 12  

        run: 
                bwa="/home/groups/CEDAR/tools/bwa" 
                ref_fa="/home/groups/CEDAR/anurpa/genomes/GRCh38.primary_assembly.genome.fa" 

                shell (""" 
                        {bwa}  mem -T 12 {ref_fa} {input[0]}  > {output}

                        """)

rule ciri2:
        input:"data/sam/bwa_map_{sample}.sam"
        output:"samples/{sample}_ciri2/{sample}_ciri2.txt"

        params:name="ciri2_{sample}", partition="long_jobs", mem="6000"
        threads:4

        run:
                ciri2="/home/groups/CEDAR/tools/CIRI2.pl"
                ref_fa="/home/groups/CEDAR/anurpa/genomes/GRCh38.primary_assembly.genome.fa"
                ref_anno="/home/groups/CEDAR/spilioto/tools/hg38_ref_all.gtf"  

                shell ("""
                        perl {ciri2} -T 12 -I {input} -O {output} -F {ref_fa} -A {ref_anno}

                        """) 

rule circ_explorer_parse:
        input: "samples/star/{sample}_bam/Chimeric.out.junction"
        output: "samples/circ_explorer/{sample}/{sample}_backspliced.junction.bed"
        params: name="circ_parse_{sample}", partition="exacloud", mem="64000"
        threads: 12
        conda:
                "../envs/circexplorer.yaml"
        shell: """CIRCexplorer2 parse -t STAR {input} -b {output}"""

rule circ_explorer_annotate:
        input: "samples/circ_explorer/{sample}/{sample}_backspliced.junction.bed"
        output: "samples/circ_explorer/{sample}/{sample}_known_circ.txt"
        params: 
                ref_fa = config["genome"],
                ref_gtf = config["ref_gtf"],
                partition="exacloud",
                mem="64000"
        threads: 4
        conda:
                "../envs/circexplorer.yaml"
        shell: """
                CIRCexplorer2 annotate -r {params.ref_gtf} -g {params.ref_fa} -b {input} -o {output}
                """
