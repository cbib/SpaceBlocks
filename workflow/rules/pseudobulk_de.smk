rule pseudobulk_de:
    """
    Pseudobulk differential expression on the concatenated multi-sample
    dataset.  Each sample is one biological replicate.

    Wildcards:
    - annot_type:     tsv_annotation | refined_annotation
    - analysis_level: by_region | by_celltype_region (| by_niche_region)

    For each level:
    - Holistic LRT: does the factor (region) explain variance?
    - Pairwise Wald: all region pairs, BH-corrected

    Input is the concatenated h5ad from integrate_samples.
    """
    input:
        adata=f"{OUTDIR_PP}/integrated_samples/concatenated.h5ad",
    output:
        results_dir=directory(
            f"{OUTDIR_PP}/pseudobulk_de/{{annot_type}}/{{analysis_level}}"
        ),
    wildcard_constraints:
        annot_type="tsv_annotation|refined_annotation|ingest_annotation",
        analysis_level="by_region|by_celltype_region|by_niche_region",
    params:
        annot_type=lambda wc: wc.annot_type,
        analysis_level=lambda wc: wc.analysis_level,
        min_cells_per_pseudobulk=ANALYSIS.get("min_cells_per_pseudobulk", 10),
        min_replicates=ANALYSIS.get("min_replicates", 3),
        de_n_genes=ANALYSIS.get("de_n_genes", 25),
    log:
        out=f"{LOGDIR}/pseudobulk_de/{{annot_type}}_{{analysis_level}}.out",
        err=f"{LOGDIR}/pseudobulk_de/{{annot_type}}_{{analysis_level}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_de/{{annot_type}}_{{analysis_level}}.tsv"
    conda:
        "../envs/pseudobulk_de_env.yaml"
    threads:
        get_resource("pseudobulk_de", "threads")
    resources:
        mem_mb=get_resource("pseudobulk_de", "mem_mb"),
        runtime=get_resource("pseudobulk_de", "runtime"),
    script:
        "../scripts/pseudobulk_de.py"
