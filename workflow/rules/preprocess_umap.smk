# workflow/rules/preprocess_umap.smk
# ─────────────────────────────────────────────────────────────────────────────
# CORE preprocessing. Now reads the UNFILTERED contract h5ad
# (config["contract"]["unfiltered_h5ad"], produced by prepare_input) instead of
# the Space Ranger tree, so it is technology-agnostic: QC filter → normalise →
# PCA → UMAP → multi-resolution Leiden, with optional precomputed-cluster reload.
# (The qupath_image output moved to rule generate_qupath_image; SR/image/geojson
# reading moved to rule prepare_input.)
# ─────────────────────────────────────────────────────────────────────────────

rule preprocess_umap:
    """Core QC / normalise / PCA / UMAP / multi-resolution Leiden on the contract h5ad."""
    input:
        unpack(_preprocess_inputs),
    output:
        adata=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
        metadata=f"{SAMPLES_DIR}/{{sample}}/metadata_{{sample}}.tsv",
        report=f"{SAMPLES_DIR}/{{sample}}/{{sample}}_report.tsv",
    params:
        sample_id=lambda wc: wc.sample,
        sample_meta=lambda wc: core_sample_meta(wc.sample),   # design columns → obs + report
        # analysis.* filtering defaults; the script may override any of these
        # per sample from the optional thresholds_tsv (absent → these defaults).
        min_counts=ANALYSIS.get("min_counts", 1),
        min_cells=ANALYSIS.get("min_cells", 3),
        min_genes=ANALYSIS.get("min_genes", 100),
        max_counts=ANALYSIS.get("max_counts", None),   # optional upper bound; None = off
        max_pct_mt=ANALYSIS.get("max_pct_mt", None),   # optional upper bound; None = off
        n_neighbors=ANALYSIS.get("n_neighbors", 10),
        n_pcs=ANALYSIS.get("n_pcs", 30),
        resolution_scan_min=ANALYSIS.get("resolution_scan_min", 0.2),
        resolution_scan_max=ANALYSIS.get("resolution_scan_max", 0.8),
        resolution_scan_step=ANALYSIS.get("resolution_scan_step", 0.1),
        random_seed=RANDOM_SEED,
        use_precomputed=USE_PRECOMPUTED,
        precomputed_metadata_dir=config.get("precomputed_metadata_dir", ""),
        region_colors=ANALYSIS.get("region_colors", {}),
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
        mem_mb=mem_mb_attempt("preprocess_umap"),
        runtime=get_resource("preprocess_umap", "runtime"),
    script:
        "../scripts/preprocess_umap.py"
