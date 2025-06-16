contrast = get_contrast

rule deseq2_init:
    input:
        counts="{project_id}_QC/{project_id}_subset_counts_pass.txt".format(project_id=config["project_id"])
    output:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        list_rds="results/diffexp/pairwise/{contrast}_all_listnorm.rds",
        ERCC_rds="results/diffexp/pairwise/{contrast}_all_ERCCnorm.rds",
        rld_out="results/diffexp/pairwise/{contrast}_rlog_dds.rds",
        list_rld_out="results/diffexp/pairwise/{contrast}_rlog_dds_listnorm.rds",
        ERCC_rld_out="results/diffexp/pairwise/{contrast}_rlog_dds_ERCCnorm.rds"
    params:
        samples=config["omic_meta_data"],
        sample_id=config["sample_id"],
        linear_model=config["linear_model"],
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule deseq2_pairwise:
    input:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        list_rds="results/diffexp/pairwise/{contrast}_all_listnorm.rds",
        ERCC_rds="results/diffexp/pairwise/{contrast}_all_ERCCnorm.rds",
        rld_out="results/diffexp/pairwise/{contrast}_rlog_dds.rds",
        list_rld_out="results/diffexp/pairwise/{contrast}_rlog_dds_listnorm.rds",
        ERCC_rld_out="results/diffexp/pairwise/{contrast}_rlog_dds_ERCCnorm.rds"
    output:
        table="results/diffexp/pairwise/{contrast}.diffexp.tsv",
        geneID_table="results/diffexp/pairwise/{contrast}.diffexp.geneID.tsv",
        ma_plot="results/diffexp/pairwise/{contrast}.ma_plot.pdf",
        p_hist="results/diffexp/pairwise/{contrast}.phist_plot.pdf",
        heatmap_plot="results/diffexp/pairwise/{contrast}.heatmap_plot.pdf",
        panel_ma="results/diffexp/pairwise/{contrast}.panel_ma.pdf",
        var_heat="results/diffexp/pairwise/{contrast}.variance_heatmap.pdf",
        pca_plot="results/diffexp/pairwise/{contrast}.pca_plot.pdf",
        list_table="results/diffexp/pairwise/{contrast}.diffexp.listnorm.tsv",
        list_geneID_table="results/diffexp/pairwise/{contrast}.diffexp.geneID.listnorm.tsv",
        list_ma_plot="results/diffexp/pairwise/{contrast}.ma_plot.listnorm.pdf",
        list_p_hist="results/diffexp/pairwise/{contrast}.phist_plot.listnorm.pdf",
        list_heatmap_plot="results/diffexp/pairwise/{contrast}.heatmap_plot.listnorm.pdf",
        list_panel_ma="results/diffexp/pairwise/{contrast}.panel_ma.listnorm.pdf",
        list_var_heat="results/diffexp/pairwise/{contrast}.variance_heatmap.listnorm.pdf",
        list_pca_plot="results/diffexp/pairwise/{contrast}.pca_plot.listnorm.pdf",
        ERCC_table="results/diffexp/pairwise/{contrast}.diffexp.ERCCnorm.tsv",
        ERCC_geneID_table="results/diffexp/pairwise/{contrast}.diffexp.geneID.ERCCnorm.tsv",
        ERCC_ma_plot="results/diffexp/pairwise/{contrast}.ma_plot.ERCCnorm.pdf",
        ERCC_p_hist="results/diffexp/pairwise/{contrast}.phist_plot.ERCCnorm.pdf",
        ERCC_heatmap_plot="results/diffexp/pairwise/{contrast}.heatmap_plot.ERCCnorm.pdf",
        ERCC_panel_ma="results/diffexp/pairwise/{contrast}.panel_ma.ERCCnorm.pdf",
        ERCC_var_heat="results/diffexp/pairwise/{contrast}.variance_heatmap.ERCCnorm.pdf",
        ERCC_pca_plot="results/diffexp/pairwise/{contrast}.pca_plot.ERCCnorm.pdf"
    params:
        contrast=get_contrast,
        linear_model=config["linear_model"],
        pca_labels=config["pca"]["labels"],
        sample_id=config["sample_id"]
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2_pairwise.R"


rule deseq2_group:
    input:
        counts="{project_id}_QC/{project_id}_subset_counts_pass.txt".format(project_id=config["project_id"])
    output:
        pca="results/diffexp/group/LRT_pca.pdf",
        sd_mean_plot="results/diffexp/group/LRT_sd_mean_plot.pdf",
        distance_plot="results/diffexp/group/LRT_distance_plot.pdf",
        heatmap_plot="results/diffexp/group/LRT_heatmap_plot.pdf",
        list_pca="results/diffexp/group/LRT_pca_listnorm.pdf",
        list_sd_mean_plot="results/diffexp/group/LRT_sd_mean_plot_listnorm.pdf",
        list_distance_plot="results/diffexp/group/LRT_distance_plot_listnorm.pdf",
        list_heatmap_plot="results/diffexp/group/LRT_heatmap_plot_listnorm.pdf",
        ERCC_pca="results/diffexp/group/LRT_pca_ERCCnorm.pdf",
        ERCC_sd_mean_plot="results/diffexp/group/LRT_sd_mean_plot_ERCCnorm.pdf",
        ERCC_distance_plot="results/diffexp/group/LRT_distance_plot_ERCC.pdf",
        ERCC_heatmap_plot="results/diffexp/group/LRT_heatmap_plot_ERCC.pdf",
        rds="results/diffexp/group/LRT_all.rds",
        rld_out="results/diffexp/group/LRT_rlog_dds.rds",
        list_rds="results/diffexp/group/LRT_all_listnorm.rds",
        list_rld_out="results/diffexp/group/LRT_rlog_dds_listnorm.rds",
        ERCC_rds="results/diffexp/group/LRT_all_ERCCnorm.rds",
        ERCC_rld_out="results/diffexp/group/LRT_rlog_dds_ERCCnorm.rds"
    params:
        pca_labels=config["pca"]["labels"],
        samples=config["omic_meta_data"],
        sample_id=config["sample_id"],
        linear_model=config["linear_model"],
        LRT=config["diffexp"]["LRT"],
        colors=config['colors']['rcolorbrewer'],
        discrete=config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_group.R"

rule deseq2_QC:
    input:
        rds="results/diffexp/group/LRT_all.rds",
        rld_out="results/diffexp/group/LRT_rlog_dds.rds",
        list_rds="results/diffexp/group/LRT_all_listnorm.rds",
        list_rld_out="results/diffexp/group/LRT_rlog_dds_listnorm.rds",
        ERCC_rds="results/diffexp/group/LRT_all_ERCCnorm.rds",
        ERCC_rld_out="results/diffexp/group/LRT_rlog_dds_ERCCnorm.rds"
    output:
        mds_plot="results/diffexp/group/MDS_plot.pdf",
        mds_table="results/diffexp/group/MDS_table.txt",
        heatmap_plot="results/diffexp/group/Heatmap_all_genes.pdf",
        sd_plot="results/diffexp/group/stdev_plot.pdf",
        rlogCounts_plot="results/diffexp/group/rlog_counts_violinPlot.pdf",
        rlogCounts_fac_plot="results/diffexp/group/rlog_counts_faceted_violinPlot.pdf",
        counts_plot="results/diffexp/group/counts_violinPlot.pdf",
        counts_fac_plot="results/diffexp/group/counts_faceted_violinPlot.pdf",
        list_mds_plot="results/diffexp/group/MDS_plot_listnorm.pdf",
        list_mds_table="results/diffexp/group/MDS_table_listnorm.txt",
        list_heatmap_plot="results/diffexp/group/Heatmap_all_genes_listnorm.pdf",
        list_sd_plot="results/diffexp/group/stdev_plot_listnorm.pdf",
        list_rlogCounts_plot="results/diffexp/group/rlog_counts_violinPlot_listnorm.pdf",
        list_rlogCounts_fac_plot="results/diffexp/group/rlog_counts_faceted_violinPlot_listnorm.pdf",
        list_counts_plot="results/diffexp/group/counts_violinPlot_listnorm.pdf",
        list_counts_fac_plot="results/diffexp/group/counts_faceted_violinPlot_listnorm.pdf",
        ERCC_mds_plot="results/diffexp/group/MDS_plot_ERCCnorm.pdf",
        ERCC_mds_table="results/diffexp/group/MDS_table_ERCCnorm.txt",
        ERCC_heatmap_plot="results/diffexp/group/Heatmap_all_genes_ERCCnorm.pdf",
        ERCC_sd_plot="results/diffexp/group/stdev_plot_ERCCnorm.pdf",
        ERCC_rlogCounts_plot="results/diffexp/group/rlog_counts_violinPlot_ERCCnorm.pdf",
        ERCC_rlogCounts_fac_plot="results/diffexp/group/rlog_counts_faceted_violinPlot_ERCCnorm.pdf",
        ERCC_counts_plot="results/diffexp/group/counts_violinPlot_ERCCnorm.pdf",
        ERCC_counts_fac_plot="results/diffexp/group/counts_faceted_violinPlot_ERCCnorm.pdf"
    params:
        sample_id=config["sample_id"],
        linear_model=config["linear_model"],
        colors=config['colors']['rcolorbrewer'],
        discrete=config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/QC.R"

rule deseq2_qplot:
    input:
        stats_table="results/diffexp/pairwise/{contrast}.diffexp.tsv",
        list_stats_table="results/diffexp/pairwise/{contrast}.diffexp.listnorm.tsv",
        ERCC_stats_table="results/diffexp/pairwise/{contrast}.diffexp.ERCCnorm.tsv"
    output:
        qplot="results/diffexp/pairwise/{contrast}.qplot.pdf",
        qhist="results/diffexp/pairwise/{contrast}.qhist.pdf",
        table="results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv",
        list_qplot="results/diffexp/pairwise/{contrast}.qplot.listnorm.pdf",
        list_qhist="results/diffexp/pairwise/{contrast}.qhist.listnorm.pdf",
        list_table="results/diffexp/pairwise/{contrast}.qvalue_diffexp.listnorm.tsv",
        ERCC_qplot="results/diffexp/pairwise/{contrast}.qplot.ERCCnorm.pdf",
        ERCC_qhist="results/diffexp/pairwise/{contrast}.qhist.ERCCnorm.pdf",
        ERCC_table="results/diffexp/pairwise/{contrast}.qvalue_diffexp.ERCCnorm.tsv"
    params:
        contrast=get_contrast,
    conda:
        "../envs/omic_qc_wf.yaml"
    script:
        "../scripts/qplot.R"


rule deseq2_density:
    input:
        rld="results/diffexp/group/LRT_rlog_dds.rds",
        list_rld="results/diffexp/group/LRT_rlog_dds_listnorm.rds",
        ERCC_rld="results/diffexp/group/LRT_rlog_dds_ERCCnorm.rds"
    output:
        density="results/diffexp/group/LRT_density_plot.pdf",
        list_density="results/diffexp/group/LRT_density_plot_listnorm.pdf",
        ERCC_density="results/diffexp/group/LRT_density_plot_ERCCnorm.pdf"
    params:
        linear_model=config["linear_model"],
        project_id=config["project_id"],
        colors=config['colors']['rcolorbrewer'],
        discrete=config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/density_plot.R"


rule GO:
    input:
        degFile="results/diffexp/pairwise/{contrast}.diffexp.tsv",
        list_degFile="results/diffexp/pairwise/{contrast}.diffexp.listnorm.tsv",
        ERCC_degFile="results/diffexp/pairwise/{contrast}.diffexp.ERCCnorm.tsv"
    output:
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"]),
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"]),
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO_listnorm.txt".format(FC = config["FC"],adjp=config["adjp"]),
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO_listnorm.txt".format(FC = config["FC"],adjp=config["adjp"]),
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO_ERCCnorm.txt".format(FC = config["FC"],adjp=config["adjp"]),
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO_ERCCnorm.txt".format(FC = config["FC"],adjp=config["adjp"])
    params:
        contrast=get_contrast,
        assembly=config["assembly"],
        printTree=config["printTree"],
        FC=config["FC"],
        adjp=config["adjp"]
    conda:
        "../envs/runGO.yaml"
    script:
        "../scripts/runGOforDESeq2.R"


rule volcano:
    input:
        degFile="results/diffexp/pairwise/{contrast}.diffexp.tsv",
        list_degFile="results/diffexp/pairwise/{contrast}.diffexp.listnorm.tsv",
        ERCC_degFile="results/diffexp/pairwise/{contrast}.diffexp.ERCCnorm.tsv"
    output:
        volcano_plot="results/diffexp/pairwise/{{contrast}}.diffexp.{adjp}.VolcanoPlot.pdf".format(adjp=config["adjp"]),
        list_volcano_plot="results/diffexp/pairwise/{{contrast}}.diffexp.{adjp}.VolcanoPlot.listnorm.pdf".format(adjp=config["adjp"]),
        ERCC_volcano_plot="results/diffexp/pairwise/{{contrast}}.diffexp.{adjp}.VolcanoPlot.ERCCnorm.pdf".format(adjp=config["adjp"])
    params:
        contrast=get_contrast,
        FC=config["FC"],
        adjp=config["adjp"]
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/RNAseq_makeVolcano.R"


rule permutation:
    input:
        counts="{project_id}_QC/{project_id}_subset_counts_pass.txt".format(project_id=config["project_id"])
    output:
        numGenes="results/diffexp/pairwise/permutationTest/{contrast}.number.diff.genes.csv",
        permList="results/diffexp/pairwise/permutationTest/{contrast}.permutation.list.csv",
        histogram="results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf",
        list_numGenes="results/diffexp/pairwise/permutationTest/{contrast}.number.diff.genes.listnorm.csv",
        list_permList="results/diffexp/pairwise/permutationTest/{contrast}.permutation.list.listnorm.csv",
        list_histogram="results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.listnorm.pdf",
        ERCC_numGenes="results/diffexp/pairwise/permutationTest/{contrast}.number.diff.genes.ERCCnorm.csv",
        ERCC_permList="results/diffexp/pairwise/permutationTest/{contrast}.permutation.list.ERCCnorm.csv",
        ERCC_histogram="results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.ERCCnorm.pdf"
    params:
        contrast=get_contrast,
        samples=config["omic_meta_data"],
        sample_id=config["sample_id"],
        linear_model=config["linear_model"]
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/permutation_test.R"


rule run_glimma:
    input:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        list_rds="results/diffexp/pairwise/{contrast}_all_listnorm.rds",
        ERCC_rds="results/diffexp/pairwise/{contrast}_all_ERCCnorm.rds"
    output:
        ma_plot="results/diffexp/glimma-plots/{contrast}.ma_plot.html",
        volcano_plot="results/diffexp/glimma-plots/{contrast}.volcano_plot.html",
        list_ma_plot="results/diffexp/glimma-plots/{contrast}.ma_plot_listnorm.html",
        list_volcano_plot="results/diffexp/glimma-plots/{contrast}.volcano_plot_listnorm.html",
        ERCC_ma_plot="results/diffexp/glimma-plots/{contrast}.ma_plot_ERCCnorm.html",
        ERCC_volcano_plot="results/diffexp/glimma-plots/{contrast}.volcano_plot_ERCCnorm.html"
    params:
        contrast=get_contrast,
        condition=config["linear_model"]
    conda:
        "../envs/glimma_env.yaml"
    script:
        "../scripts/run_glimma.R"


rule run_glimma_mds:
    input:
        rds="results/diffexp/group/LRT_all.rds",
        list_rds="results/diffexp/group/LRT_all_listnorm.rds",
        ERCC_rds="results/diffexp/group/LRT_all_ERCCnorm.rds"
    output:
        mds_plot="results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
        list_mds_plot="results/diffexp/glimma-plots/{project_id}.mds_plot_listnorm.html".format(project_id=project_id),
        ERCC_mds_plot="results/diffexp/glimma-plots/{project_id}.mds_plot_ERCCnorm.html".format(project_id=project_id)
    params:
        project_id=config["project_id"],
    conda:
        "../envs/glimma_env.yaml"
    script:
        "../scripts/run_glimma_mds.R"
