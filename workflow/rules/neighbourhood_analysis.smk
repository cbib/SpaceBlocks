rule neighbourhood_analysis:
    """Spatial neighbourhood analysis per sample × annot_type."""
    input:
        adata=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_annotated.h5ad",
    output:
        results_dir=directory(f"{SAMPLES_DIR}/{{sample}}/neighbourhood_analysis/{{annot_type}}"),
    params:
        sample_id=lambda wc: wc.sample,
        annot_type=lambda wc: wc.annot_type,
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
