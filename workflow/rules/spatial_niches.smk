rule spatial_niches:
    """
    Cross-sample spatial niche / domain identification with BANKSY.

    Runs jointly on all preprocessed samples (staggered coords → per-sample
    spatial graphs → BANKSY matrix → PCA → Harmony across samples → Leiden),
    producing one spatial_niche label per cell.  The per-sample niche TSVs are
    consumed by annotate_cells, which injects the `spatial_niche` obs column so
    every downstream rule picks it up automatically (find_niche_column).

    Runs before annotation, independently of integrate_samples, so the niche
    method can be swapped (or precomputed externally) without touching the rest
    of the pipeline.
    """
    input:
        adatas=expand(f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
                      sample=SAMPLE_IDS),
    output:
        concatenated=f"{OUTDIR_PP}/spatial_niches/spatial_niches_concatenated.h5ad",
        niche_tsvs=expand(f"{OUTDIR_PP}/spatial_niches/tsv/niche_{{sample}}.tsv",
                          sample=SAMPLE_IDS),
        plots_dir=directory(f"{OUTDIR_PP}/spatial_niches/plots"),
    params:
        sample_ids=SAMPLE_IDS,
        random_seed=RANDOM_SEED,
        banksy_lambda=SPATIAL_NICHES.get("lambda", 0.8),
        num_neighbours=SPATIAL_NICHES.get("num_neighbours", 18),
        max_m=SPATIAL_NICHES.get("max_m", 1),
        nbr_weight_decay=SPATIAL_NICHES.get("nbr_weight_decay", "scaled_gaussian"),
        use_hvg=SPATIAL_NICHES.get("use_hvg", False),
        n_top_genes=SPATIAL_NICHES.get("n_top_genes", ANALYSIS.get("n_top_genes", 2000)),
        n_pcs=SPATIAL_NICHES.get("n_pcs", ANALYSIS.get("n_pcs", 30)),
        niche_resolution=SPATIAL_NICHES.get("niche_resolution", 0.3),
        resolution_scan_min=SPATIAL_NICHES.get("resolution_scan_min", 0.2),
        resolution_scan_max=SPATIAL_NICHES.get("resolution_scan_max", 1.0),
        resolution_scan_step=SPATIAL_NICHES.get("resolution_scan_step", 0.2),
        run_resolution_scan=SPATIAL_NICHES.get("run_resolution_scan", False),
        compute_umap=SPATIAL_NICHES.get("compute_umap", False),
        cluster_n_neighbours=SPATIAL_NICHES.get("cluster_n_neighbours", 15),
        n_iterations=SPATIAL_NICHES.get("n_iterations", 2),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        dpi=ANALYSIS.get("plot_dpi", 300),
    log:
        out=f"{LOGDIR}/spatial_niches/spatial_niches.out",
        err=f"{LOGDIR}/spatial_niches/spatial_niches.err",
    benchmark:
        f"{LOGDIR}/benchmarks/spatial_niches/spatial_niches.tsv"
    conda:
        "../envs/spatial_niches.yaml"
    threads:
        get_resource("spatial_niches", "threads")
    resources:
        mem_mb=get_resource("spatial_niches", "mem_mb"),
        runtime=get_resource("spatial_niches", "runtime"),
    script:
        "../scripts/spatial_niches.py"
