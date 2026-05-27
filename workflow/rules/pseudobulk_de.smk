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
