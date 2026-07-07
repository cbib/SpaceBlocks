rule neighbourhood_analysis:
    """
    Spatial neighbourhood analysis, run for both annotation types.

    The {annot_type} wildcard is either 'tsv_annotation' or 'refined_annotation'.
    """
    input:
        adata=rules.annotate_cells.output.adata_annot,
    output:
        results_dir=directory(f"{SAMPLES_DIR}/{{sample}}/neighbourhood_analysis/{{annot_type}}"),
    params:
        sample_id=lambda wc: wc.sample,
        annot_type=lambda wc: wc.annot_type,
        annotation_colors=config.get("annotation_colors", {}),
    log:
        out=f"{LOGDIR}/neighbourhood_analysis/{{sample}}_{{annot_type}}.out",
        err=f"{LOGDIR}/neighbourhood_analysis/{{sample}}_{{annot_type}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/neighbourhood_analysis/{{sample}}_{{annot_type}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("neighbourhood_analysis", "threads")
    resources:
        mem_mb=get_resource("neighbourhood_analysis", "mem_mb"),
        runtime=get_resource("neighbourhood_analysis", "runtime"),
    script:
        "../scripts/neighbourhood_analysis.py"
