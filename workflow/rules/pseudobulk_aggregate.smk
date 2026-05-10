rule pseudobulk_aggregate:
    """
    Aggregate single cells into pseudobulk count matrices using decoupler.

    Produces TSV count matrices and metadata per analysis level,
    ready to be consumed by the R-based pseudobulk_de rule.
    """
    input:
        adata=f"{OUTDIR_PP}/integrated_samples/concatenated.h5ad",
    output:
        agg_dir=directory(
            f"{OUTDIR_PP}/pseudobulk/{{annot_type}}/{{analysis_level}}/matrices"
        ),
    wildcard_constraints:
        annot_type="tsv_annotation|refined_annotation|ingest_annotation",
        analysis_level="by_region|by_celltype_region|by_niche_region",
    params:
        annot_type=lambda wc: wc.annot_type,
        analysis_level=lambda wc: wc.analysis_level,
        min_cells_per_pseudobulk=ANALYSIS.get("min_cells_per_pseudobulk", 10),
        min_counts_per_pseudobulk=ANALYSIS.get("min_counts_per_pseudobulk", 1000),
    log:
        out=f"{LOGDIR}/pseudobulk_aggregate/{{annot_type}}_{{analysis_level}}.out",
        err=f"{LOGDIR}/pseudobulk_aggregate/{{annot_type}}_{{analysis_level}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_aggregate/{{annot_type}}_{{analysis_level}}.tsv"
    conda:
        "../envs/pseudobulk_aggregate.yaml"
    threads:
        get_resource("pseudobulk_aggregate", "threads")
    resources:
        mem_mb=get_resource("pseudobulk_aggregate", "mem_mb"),
        runtime=get_resource("pseudobulk_aggregate", "runtime"),
    script:
        "../scripts/pseudobulk_aggregate.py"
