rule sample_report:
    """
    Generate a multi-page PDF report with separate pages per sample.

    Each sample contains UMAP/spatial panels for every available annotation,
    plus niche/region summaries, primary-annotation composition, and a marker
    dotplot. External annotation is prioritised over ingest when configured.
    """
    input:
        annotated=expand(rules.annotate_cells.output.adata_annot, sample=SAMPLE_IDS),
    output:
        report=f"{OUTDIR_PP}/integrated_samples/samples_report.pdf",
    log:
        out=f"{LOGDIR}/sample_report/sample_report.out",
        err=f"{LOGDIR}/sample_report/sample_report.err",
    benchmark:
        f"{LOGDIR}/benchmarks/sample_report/sample_report.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads: get_resource("sample_report", "threads")
    resources:
        mem_mb=mem_mb_attempt("sample_report"),
        runtime=get_resource("sample_report", "runtime"),
    params:
        sample_ids=SAMPLE_IDS,
        primary_annotation_column=DEFAULT_ANNOT_COL,
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        dpi=ANALYSIS.get("plot_dpi", 300),
        umap_point_size=ANALYSIS.get("umap_point_size", 2),
        spatial_point_size=ANALYSIS.get("spatial_point_size", 20),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
    script:
        "../scripts/sample_report.py"
