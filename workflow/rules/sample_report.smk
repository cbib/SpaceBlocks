rule sample_report:
    """
    Generate a multi-page PDF report with one page per sample.

    Each page contains: UMAP (clusters + annotation), spatial plot
    (clusters + annotation), dotplot with cell type markers, and
    a barplot with cell proportions per region.
    """
    input:
        annotated=expand(
            f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_annotated.h5ad",
            sample=SAMPLE_IDS,
        ),
        annotation_markers=config.get("annotation_markers", ""),
    output:
        report=f"{OUTDIR_PP}/integrated_samples/samples_report.pdf",
    params:
        sample_ids=SAMPLE_IDS,
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
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
        mem_mb=get_resource("sample_report", "mem_mb"),
        runtime=get_resource("sample_report", "runtime"),
    script:
        "../scripts/sample_report.py"
