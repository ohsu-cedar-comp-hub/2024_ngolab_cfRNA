rule QC_init:
    input:
        counts="data/{project_id}_counts.txt".format(project_id=config["project_id"]),
        coverage="results/tables/read_coverage.txt",
        st_stats="results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id=config["project_id"])
    output:
        exon_fraction="{project_id}_QC/Exon.Fraction.pdf".format(project_id=config["project_id"]),
        input_reads="{project_id}_QC/Input.Reads.Barplot.pdf".format(project_id=config["project_id"]),
        unique_reads="{project_id}_QC/Unique.Reads.Barplot.pdf".format(project_id=config["project_id"]),
        lib_size="{project_id}_QC/Scatter.Plot.LibSize.ExonFraction.pdf".format(project_id=config["project_id"]),
        md_stats="{project_id}_QC/{project_id}_Metadata_W_QC_Metrics.txt".format(project_id=config["project_id"]),
        batch_pca="{project_id}_QC/{project_id}_PCA_potential_batch_grid.pdf".format(project_id=config["project_id"]),
        sub_counts="{project_id}_QC/{project_id}_subset_counts_RPM.txt".format(project_id=config["project_id"]),
        subset_cts="{project_id}_QC/{project_id}_subset_counts_pass.txt".format(project_id=config["project_id"])
    params:
        samples=config["omic_meta_data"],
        sample_id = config["sample_id"],
        ens2geneID=config["ens2geneID"]
    conda:
        "../envs/QC.yaml"
    script:
        "../scripts/QC_init.R"

