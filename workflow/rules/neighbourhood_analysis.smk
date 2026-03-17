rule neighbourhood_analysis:
    """
    Spatial neighbourhood analysis using squidpy.

    Computes neighbourhood enrichment and co-occurrence scores based on
    cell-type annotations.  Requires the same cluster annotation TSV
    as the pseudobulk_de rule.
    """
    input:
        adata=f"{OUTDIR_PP}/{{sample}}/adata_{{sample}}.h5ad",
        cluster_annotations=config.get("cluster_annotations", ""),
    output:
        results_dir=directory(f"{OUTDIR_PP}/{{sample}}/neighbourhood_analysis"),
    params:
        sample_id=lambda wc: wc.sample,
    log:
        out=f"{LOGDIR}/neighbourhood_analysis/{{sample}}.out",
        err=f"{LOGDIR}/neighbourhood_analysis/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/neighbourhood_analysis/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("neighbourhood_analysis", "threads")
    resources:
        mem_mb=get_resource("neighbourhood_analysis", "mem_mb"),
        runtime=get_resource("neighbourhood_analysis", "runtime"),
    script:
        "../scripts/neighbourhood_analysis.py"
