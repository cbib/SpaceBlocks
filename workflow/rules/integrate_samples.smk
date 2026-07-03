rule integrate_samples:
    """
    Integrate all annotated samples into a multi-sample object.

    Produces three h5ad files in integrated_samples/:
    - concatenated.h5ad:  simple concatenation (no batch correction)
    - harmony_integrated.h5ad:  Harmony-corrected PCA + UMAP
    - sketched.h5ad:  geosketched 25% subset with projected clusters
    """
    input:
        annotated=expand(
            f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_annotated.h5ad",
            sample=SAMPLE_IDS,
        ),
    output:
        concatenated=f"{OUTDIR_PP}/integrated_samples/concatenated.h5ad",
        harmony=f"{OUTDIR_PP}/integrated_samples/harmony_integrated.h5ad",
        sketched=f"{OUTDIR_PP}/integrated_samples/sketched.h5ad",
    params:
        sample_ids=SAMPLE_IDS,
        n_neighbors=ANALYSIS.get("n_neighbors", 10),
        sketch_fraction=ANALYSIS.get("sketch_fraction", 0.25),
        random_seed=RANDOM_SEED,
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        dpi=ANALYSIS.get("plot_dpi", 300),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
        extra_annot_columns=EXTRA_ANNOT_COLUMNS,
        sample_colors=SAMPLE_COLORS,
        integrate_key=INTEGRATE_KEY,
    log:
        out=f"{LOGDIR}/integrate_samples/integrate.out",
        err=f"{LOGDIR}/integrate_samples/integrate.err",
    benchmark:
        f"{LOGDIR}/benchmarks/integrate_samples/integrate.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("integrate_samples", "threads")
    resources:
        mem_mb=get_resource("integrate_samples", "mem_mb"),
        runtime=get_resource("integrate_samples", "runtime"),
    script:
        "../scripts/integrate_samples.py"
