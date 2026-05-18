rule leiden_analysis:
    """
    Per-resolution visualisation.

    Reads the adata (which already contains leiden_{res} in obs), selects
    the appropriate column, and produces all plots.  No h5ad is saved
    per resolution — the single adata from preprocess_umap has everything.
    """
    input:
        adata=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
        cell_markers=config["snakemake_cell_markers"],
    output:
        res_dir=directory(f"{SAMPLES_DIR}/{{sample}}/leiden_resolution_{{resolution}}"),
    params:
        sample_id=lambda wc: wc.sample,
        resolution=lambda wc: wc.resolution,
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
        resolution_scan_min=ANALYSIS.get("resolution_scan_min", 0.2),
        resolution_scan_max=ANALYSIS.get("resolution_scan_max", 0.8),
        resolution_scan_step=ANALYSIS.get("resolution_scan_step", 0.1),
    log:
        out=f"{LOGDIR}/leiden_analysis/{{sample}}_res{{resolution}}.out",
        err=f"{LOGDIR}/leiden_analysis/{{sample}}_res{{resolution}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/leiden_analysis/{{sample}}_res{{resolution}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("leiden_analysis", "threads")
    resources:
        mem_mb=get_resource("leiden_analysis", "mem_mb"),
        runtime=get_resource("leiden_analysis", "runtime"),
    script:
        "../scripts/leiden_analysis.py"
