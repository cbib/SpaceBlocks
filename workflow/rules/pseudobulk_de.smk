rule pseudobulk_de:
    """Explicit pairwise Wald contrasts and optional omnibus LRT with R/DESeq2."""
    input:
        agg_dir=rules.pseudobulk_aggregate.output.agg_dir,
    output:
        results_dir=directory(
            f"{OUTDIR_PP}/pseudobulk/{{annot_type}}/{{analysis_name}}/de_results"
        ),
    log:
        out=f"{LOGDIR}/pseudobulk_de/{{annot_type}}_{{analysis_name}}.out",
        err=f"{LOGDIR}/pseudobulk_de/{{annot_type}}_{{analysis_name}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_de/{{annot_type}}_{{analysis_name}}.tsv"
    wildcard_constraints:
        annot_type="tsv_annotation|external_annotation|ingest_annotation|refined_annotation",
        analysis_name="[A-Za-z0-9_-]+",
    conda:
        "../envs/pseudobulk_de.yaml"
    threads: get_resource("pseudobulk_de", "threads")
    resources:
        mem_mb=mem_mb_attempt("pseudobulk_de"),
        runtime=get_resource("pseudobulk_de", "runtime"),
    params:
        analysis_name=lambda wc: wc.analysis_name,
        min_replicates=ANALYSIS.get("min_replicates", 3),
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
        padj_threshold=ANALYSIS.get("padj_threshold", 0.05),
        lfc_threshold=ANALYSIS.get("lfc_threshold", 0.5),
        dpi=ANALYSIS.get("plot_dpi", 300),
        group_by_columns=lambda wc: pseudobulk_analysis(wc.analysis_name)["group_by"],
        paired_by=lambda wc: pseudobulk_analysis(wc.analysis_name)["paired_by"],
        covariates=lambda wc: pseudobulk_analysis(wc.analysis_name)["covariates"],
        contrast_names=lambda wc: pseudobulk_contrast_param(wc.analysis_name, "names"),
        contrast_numerators=lambda wc: pseudobulk_contrast_param(
            wc.analysis_name, "numerator"
        ),
        contrast_denominators=lambda wc: pseudobulk_contrast_param(
            wc.analysis_name, "denominator"
        ),
        lrt_enabled=lambda wc: bool(
            pseudobulk_analysis(wc.analysis_name)["lrt"].get("enabled", False)
        ),
        group_color_names=lambda wc: list(
            pseudobulk_group_palette(wc.analysis_name).keys()
        ),
        group_color_values=lambda wc: list(
            pseudobulk_group_palette(wc.analysis_name).values()
        ),
        condition_level_order=lambda wc: pseudobulk_condition_order(wc.analysis_name),
        extra_annot_columns=EXTRA_ANNOT_COLUMNS,
        extra_anno_col_names=EXTRA_ANNO_COLS,
        extra_anno_values=EXTRA_ANNO_VALS,
        extra_anno_colors=EXTRA_ANNO_COLORS,
    script:
        "../scripts/pseudobulk_de.R"
