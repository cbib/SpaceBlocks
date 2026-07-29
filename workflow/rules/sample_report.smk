rule sample_report:
    """
    Generate a multi-page PDF report with separate pages per sample.

    Each page contains: (1) UMAP (clusters + annotation), spatial plot
    (clusters + annotation), (2-3) barplots with cell proportions per region, cluster and niche,
    and (4) a dotplot with the top markers per cell type.
    """
    input:
        annotated=expand(rules.annotate_cells.output.adata_annot, sample=SAMPLE_IDS),
    output:
        report=f"{OUTDIR_PP}/integrated_samples/samples_report.pdf",
    params:
        sample_ids=SAMPLE_IDS,
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        dpi=ANALYSIS.get("plot_dpi", 300),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
    log:
        out=f"{LOGDIR}/sample_report/sample_report.out",
        err=f"{LOGDIR}/sample_report/sample_report.err",
    benchmark:
        f"{LOGDIR}/benchmarks/sample_report/sample_report.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("sample_report", "threads")
    resources:
        mem_mb=mem_mb_attempt("sample_report"),
        runtime=get_resource("sample_report", "runtime"),
    script:
        "../scripts/sample_report.py"
