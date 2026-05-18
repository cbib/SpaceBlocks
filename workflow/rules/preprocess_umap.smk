rule preprocess_umap:
    """
    Load Space Ranger outputs, QC, normalise, PCA, Harmony, UMAP, all Leiden.

    Discovers files from the Space Ranger output directory rather than
    depending on hardcoded paths, since the tree varies across versions.
    """
    input:
        sr_done=f"{OUTDIR_SR}/{{sample}}/.done",
    output:
        adata=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
        metadata=f"{SAMPLES_DIR}/{{sample}}/metadata_{{sample}}.tsv",
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_tissue_hires_image.png",
    params:
        sample_id=lambda wc: wc.sample,
        sr_outdir=f"{OUTDIR_SR}/{{sample}}",
        geojson_path=GEOJ_DIR,
        min_counts=ANALYSIS.get("min_counts", 1),
        min_cells=ANALYSIS.get("min_cells", 3),
        min_genes=ANALYSIS.get("min_genes", 100),
        n_neighbors=ANALYSIS.get("n_neighbors", 10),
        n_pcs=ANALYSIS.get("n_pcs", 30),
        resolution_scan_min=ANALYSIS.get("resolution_scan_min", 0.2),
        resolution_scan_max=ANALYSIS.get("resolution_scan_max", 0.8),
        resolution_scan_step=ANALYSIS.get("resolution_scan_step", 0.1),
        random_seed=RANDOM_SEED,
        use_precomputed=USE_PRECOMPUTED,
    log:
        out=f"{LOGDIR}/preprocess_umap/{{sample}}.out",
        err=f"{LOGDIR}/preprocess_umap/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/preprocess_umap/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("preprocess_umap", "threads")
    resources:
        mem_mb=get_resource("preprocess_umap", "mem_mb"),
        runtime=get_resource("preprocess_umap", "runtime"),
    script:
        "../scripts/preprocess_umap.py"
