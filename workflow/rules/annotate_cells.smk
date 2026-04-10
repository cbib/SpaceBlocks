def _annotate_input_adata(wc):
    """
    If ingest_ref is configured, use the ingested h5ad (which carries
    cell_type_ingest in obs).  Otherwise use the base preprocessed adata.
    """
    if config.get("ingest_ref", ""):
        return f"{OUTDIR_PP}/{wc.sample}/adata_{wc.sample}_ingested.h5ad"
    return f"{OUTDIR_PP}/{wc.sample}/adata_{wc.sample}.h5ad"


rule annotate_cells:
    """
    Annotate cells from cluster annotations TSV, then refine using
    scaled marker expression from annotation_markers.tsv.

    If ingest_ref was run, the input adata already contains
    obs['cell_type_ingest'] from the reference transfer.
    """
    input:
        adata=_annotate_input_adata,
        cluster_annotations=config.get("cluster_annotations", ""),
        annotation_markers=config.get("annotation_markers", ""),
    output:
        adata_annot=f"{OUTDIR_PP}/{{sample}}/adata_{{sample}}_annotated.h5ad",
        plots_dir=directory(f"{OUTDIR_PP}/{{sample}}/annotation"),
    params:
        sample_id=lambda wc: wc.sample,
        marker_threshold=ANALYSIS.get("marker_refinement_threshold", 1.0),
        min_markers_expressed=ANALYSIS.get("min_markers_expressed", 2),
        min_cells_per_type=ANALYSIS.get("min_cells_per_type", 15),
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
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
