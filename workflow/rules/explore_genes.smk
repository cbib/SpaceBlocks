rule explore_genes_integrated:
    """
    Gene / signature exploration – integrated plots (PNGs).

    Loads the Harmony-integrated h5ad once, computes AUCell scores for
    all signatures, writes shared expression ranges, and produces PNGs
    in per-entry subdirectories under gene_exploration/{entry}/Integrated/.
    """
    input:
        integrated=f"{OUTDIR_PP}/integrated_samples/harmony_integrated.h5ad",
        queries=GENE_EXPLORATION.get("queries", ""),
    output:
        ranges=f"{OUTDIR_PP}/gene_exploration/expression_ranges.tsv",
        done=touch(f"{OUTDIR_PP}/gene_exploration/.integrated_done"),
    params:
        outdir=f"{OUTDIR_PP}/gene_exploration",
        annot_key=GENE_EXPLORATION.get("annot_key", "cell_type_tsv"),
        aucell_fraction=GENE_EXPLORATION.get("aucell_max_rank_fraction", 0.05),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
        dpi=GENE_EXPLORATION.get("dpi", ANALYSIS.get("plot_dpi", 300)),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        extra_annot_columns=EXTRA_ANNOT_COLUMNS,
        sample_colors=SAMPLE_COLORS,
    log:
        out=f"{LOGDIR}/explore_genes/integrated.out",
        err=f"{LOGDIR}/explore_genes/integrated.err",
    benchmark:
        f"{LOGDIR}/benchmarks/explore_genes/integrated.tsv"
    conda:
        "../envs/pseudobulk_aggregate.yaml"
    threads:
        get_resource("explore_genes_integrated", "threads")
    resources:
        mem_mb=get_resource("explore_genes_integrated", "mem_mb"),
        runtime=get_resource("explore_genes_integrated", "runtime"),
    script:
        "../scripts/explore_genes_integrated.py"


rule explore_genes_sample:
    """
    Gene / signature exploration – per-sample plots (PNGs).

    Processes a single sample: composite spatial PNG (expression +
    cell type + region) and composite dotplot PNGs (3 views stacked).
    Parallelised across samples via Snakemake wildcard.
    Output organised by entry: gene_exploration/{entry}/Spatial/ and
    gene_exploration/{entry}/Dotplots/.
    """
    input:
        adata=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_annotated.h5ad",
        ranges=f"{OUTDIR_PP}/gene_exploration/expression_ranges.tsv",
        queries=GENE_EXPLORATION.get("queries", ""),
    output:
        done=touch(f"{OUTDIR_PP}/gene_exploration/.sentinels/{{sample}}.done"),
    params:
        outdir=f"{OUTDIR_PP}/gene_exploration",
        sample_id=lambda wc: wc.sample,
        annot_key=GENE_EXPLORATION.get("annot_key", "cell_type_tsv"),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
        aucell_fraction=GENE_EXPLORATION.get("aucell_max_rank_fraction", 0.05),
        dpi=GENE_EXPLORATION.get("dpi", ANALYSIS.get("plot_dpi", 300)),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
    log:
        out=f"{LOGDIR}/explore_genes/{{sample}}.out",
        err=f"{LOGDIR}/explore_genes/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/explore_genes/{{sample}}.tsv"
    conda:
        "../envs/pseudobulk_aggregate.yaml"
    threads:
        get_resource("explore_genes_sample", "threads")
    resources:
        mem_mb=get_resource("explore_genes_sample", "mem_mb"),
        runtime=get_resource("explore_genes_sample", "runtime"),
    script:
        "../scripts/explore_genes_sample.py"
