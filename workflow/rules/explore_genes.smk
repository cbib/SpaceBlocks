rule explore_genes_integrated:
    """
    Gene / signature exploration – integrated plots.

    Loads the Harmony-integrated h5ad once, computes AUCell scores for
    all signatures, writes shared expression ranges, and produces one
    integrated PDF per gene/signature entry.
    """
    input:
        integrated=f"{OUTDIR_PP}/integrated_samples/harmony_integrated.h5ad",
        queries=GENE_EXPLORATION.get("queries", ""),
    output:
        ranges=f"{OUTDIR_PP}/gene_exploration/expression_ranges.tsv",
        outdir=directory(f"{OUTDIR_PP}/gene_exploration/integrated"),
    params:
        annot_key=GENE_EXPLORATION.get("annot_key", "cell_type_tsv"),
        aucell_fraction=GENE_EXPLORATION.get("aucell_max_rank_fraction", 0.05),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
        dpi=GENE_EXPLORATION.get("dpi", 300),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
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


rule explore_genes_spatial:
    """
    Gene / signature exploration – per-sample spatial plots.

    Loads each annotated sample sequentially (one at a time), computes
    per-sample AUCell scores, and appends spatial pages into one PDF
    per gene/signature entry.  Memory bounded to a single sample.
    """
    input:
        annotated=expand(
            f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_annotated.h5ad",
            sample=SAMPLE_IDS,
        ),
        ranges=f"{OUTDIR_PP}/gene_exploration/expression_ranges.tsv",
        queries=GENE_EXPLORATION.get("queries", ""),
    output:
        done=touch(f"{OUTDIR_PP}/gene_exploration/spatial/.done"),
    params:
        sample_ids=SAMPLE_IDS,
        annot_key=GENE_EXPLORATION.get("annot_key", "cell_type_tsv"),
        aucell_fraction=GENE_EXPLORATION.get("aucell_max_rank_fraction", 0.05),
        dpi=GENE_EXPLORATION.get("dpi", 300),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
    log:
        out=f"{LOGDIR}/explore_genes/spatial.out",
        err=f"{LOGDIR}/explore_genes/spatial.err",
    benchmark:
        f"{LOGDIR}/benchmarks/explore_genes/spatial.tsv"
    conda:
        "../envs/pseudobulk_aggregate.yaml"
    threads:
        get_resource("explore_genes_spatial", "threads")
    resources:
        mem_mb=get_resource("explore_genes_spatial", "mem_mb"),
        runtime=get_resource("explore_genes_spatial", "runtime"),
    script:
        "../scripts/explore_genes_spatial.py"
