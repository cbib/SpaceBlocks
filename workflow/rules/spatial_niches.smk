def _spatial_niches_inputs(wc):
    """Inputs for spatial_niches. When use_precomputed is set with an external
    niche_dir, the per-sample niche_{sample}.tsv are added as MANDATORY inputs:
    Snakemake then aborts the rule (Missing input files) if any is absent, before
    the job starts. This is complementary to the in-script safeguard, which also
    catches empty TSVs and partial cell coverage (which Snakemake cannot see).
    When niche_dir is empty the rule reloads from its own output dir, so nothing
    is added here (no circular dependency) and the script check applies alone."""
    inputs = {
        "adatas": expand(f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
                         sample=SAMPLE_IDS),
    }
    sn = config.get("spatial_niches", {})
    if sn.get("use_precomputed", False):
        d = sn.get("niche_dir", "")
        if d:
            inputs["precomputed_niches"] = [
                os.path.join(d, f"niche_{s}.tsv") for s in SAMPLE_IDS
            ]
    return inputs


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
        unpack(_spatial_niches_inputs),
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
        compute_umap=SPATIAL_NICHES.get("compute_umap", False),
        cluster_n_neighbours=SPATIAL_NICHES.get("cluster_n_neighbours", 15),
        n_iterations=SPATIAL_NICHES.get("n_iterations", 2),
        use_precomputed=SPATIAL_NICHES.get("use_precomputed", False),
        niche_dir=SPATIAL_NICHES.get("niche_dir", ""),
        refine=SPATIAL_NICHES.get("refine", False),
        refine_num_neigh=SPATIAL_NICHES.get("refine_num_neigh", 6),
        refine_iterations=SPATIAL_NICHES.get("refine_iterations", 1),
        refine_auto=SPATIAL_NICHES.get("refine_auto", False),
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
