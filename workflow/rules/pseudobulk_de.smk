rule pseudobulk_de:
    """
    DE with R DESeq2 (Wald + LRT), EnhancedVolcano, ComplexHeatmap, DEGpatterns.
    Region levels and colors from config control ordering and visualization.
    """
    input:
        agg_dir=f"{OUTDIR_PP}/pseudobulk/{{annot_type}}/{{analysis_level}}/aggregated",
    output:
        results_dir=directory(
            f"{OUTDIR_PP}/pseudobulk/{{annot_type}}/{{analysis_level}}/de_results"
        ),
    wildcard_constraints:
        annot_type="tsv_annotation|refined_annotation|ingest_annotation",
        analysis_level="by_region|by_celltype_region|by_niche_region",
    params:
        annot_type=lambda wc: wc.annot_type,
        analysis_level=lambda wc: wc.analysis_level,
        min_replicates=ANALYSIS.get("min_replicates", 3),
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
        padj_threshold=ANALYSIS.get("padj_threshold", 0.05),
        lfc_threshold=ANALYSIS.get("lfc_threshold", 0.5),
        region_levels=REGION_LEVELS,
        # Pass colors as two parallel lists for safe Python→R conversion
        region_color_names=list(ANALYSIS.get("region_colors", {}).keys()),
        region_color_values=list(ANALYSIS.get("region_colors", {}).values()),
        # Design columns to annotate heatmaps with (chosen via config design.columns);
        # palettes flattened into parallel lists (grey fallback applied in R).
        extra_annot_columns=EXTRA_ANNOT_COLUMNS,
        extra_anno_col_names=EXTRA_ANNO_COLS,
        extra_anno_values=EXTRA_ANNO_VALS,
        extra_anno_colors=EXTRA_ANNO_COLORS,
    log:
        out=f"{LOGDIR}/pseudobulk_de/{{annot_type}}_{{analysis_level}}.out",
        err=f"{LOGDIR}/pseudobulk_de/{{annot_type}}_{{analysis_level}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_de/{{annot_type}}_{{analysis_level}}.tsv"
    conda:
        "../envs/pseudobulk_de.yaml"
    threads:
        get_resource("pseudobulk_de", "threads")
    resources:
        mem_mb=get_resource("pseudobulk_de", "mem_mb"),
        runtime=get_resource("pseudobulk_de", "runtime"),
    script:
        "../scripts/pseudobulk_de.R"
