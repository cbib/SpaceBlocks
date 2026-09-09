rule pseudobulk_aggregate:
    """
    Aggregate single cells into pseudobulk count matrices using a named analysis.

    Output structure:
      {annot_type}/{analysis_name}/aggregated/
        ├── matrices/     count TSVs
        ├── metadata/     metadata TSVs
        ├── plots/        sample-level QC
        └── manifest.tsv
    """
    input:
        adata=rules.integrate_samples.output.concatenated,
    output:
        agg_dir=directory(
            f"{OUTDIR_PP}/pseudobulk/{{annot_type}}/{{analysis_name}}/aggregated"
        ),
    log:
        out=f"{LOGDIR}/pseudobulk_aggregate/{{annot_type}}_{{analysis_name}}.out",
        err=f"{LOGDIR}/pseudobulk_aggregate/{{annot_type}}_{{analysis_name}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_aggregate/{{annot_type}}_{{analysis_name}}.tsv"
    wildcard_constraints:
        annot_type="tsv_annotation|external_annotation|ingest_annotation|refined_annotation",
        analysis_name="[A-Za-z0-9_-]+",
    conda:
        "../envs/pseudobulk_aggregate.yaml"
    threads: get_resource("pseudobulk_aggregate", "threads")
    resources:
        mem_mb=mem_mb_attempt("pseudobulk_aggregate"),
        runtime=get_resource("pseudobulk_aggregate", "runtime"),
    params:
        annot_type=lambda wc: wc.annot_type,
        analysis_spec=lambda wc: pseudobulk_analysis(wc.analysis_name),
        min_cells_per_pseudobulk=ANALYSIS.get("min_cells_per_pseudobulk", 10),
        min_counts_per_pseudobulk=ANALYSIS.get("min_counts_per_pseudobulk", 1000),
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        dpi=ANALYSIS.get("plot_dpi", 300),
        extra_annot_columns=EXTRA_ANNOT_COLUMNS,
        sample_colors=SAMPLE_COLORS,
    script:
        "../scripts/pseudobulk_aggregate.py"
