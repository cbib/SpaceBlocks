def _annotate_input_adata(wc):
    if config.get("ingest_ref", ""):
        return f"{SAMPLES_DIR}/{wc.sample}/adata_{wc.sample}_ingested.h5ad"
    return f"{SAMPLES_DIR}/{wc.sample}/adata_{wc.sample}.h5ad"


def _annotate_niche_input(wc):
    """Spatial-niche TSV for this sample (barcode → spatial_niche), injected as
    the spatial_niche obs column. Always the spatial_niches rule output: the rule
    itself handles precomputed reload (use_precomputed + niche_dir), so annotate
    consumes a single canonical path. Returns [] when niche identification is
    disabled, so annotate_cells does not depend on it."""
    sn = config.get("spatial_niches", {})
    if not sn.get("enabled", False):
        return []
    return f"{OUTDIR_PP}/spatial_niches/tsv/niche_{wc.sample}.tsv"


rule annotate_cells:
    """
    Annotate cells: TSV cluster mapping + optional external annotation.
    If use_precomputed_clusters, reads leiden from metadata TSV.
    """
    input:
        adata=_annotate_input_adata,
        metadata=f"{SAMPLES_DIR}/{{sample}}/metadata_{{sample}}.tsv",
        cluster_annotations=config.get("cluster_annotations", ""),
        spatial_niche=_annotate_niche_input,
    output:
        adata_annot=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_annotated.h5ad",
        plots_dir=directory(f"{SAMPLES_DIR}/{{sample}}/annotation"),
    params:
        sample_id=lambda wc: wc.sample,
        min_cells_per_type=ANALYSIS.get("min_cells_per_type", 15),
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
        use_precomputed=USE_PRECOMPUTED,
        external_annotation=config.get("external_annotation", {}),
        precomputed_metadata_dir=config.get("precomputed_metadata_dir", ""),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        dpi=ANALYSIS.get("plot_dpi", 300),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
    log:
        out=f"{LOGDIR}/annotate_cells/{{sample}}.out",
        err=f"{LOGDIR}/annotate_cells/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/annotate_cells/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("annotate_cells", "threads")
    resources:
        mem_mb=get_resource("annotate_cells", "mem_mb"),
        runtime=get_resource("annotate_cells", "runtime"),
    script:
        "../scripts/annotate_cells.py"
