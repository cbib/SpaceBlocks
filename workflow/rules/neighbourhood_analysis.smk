rule neighbourhood_analysis:
    """
    Spatial neighbourhood analysis, run for both annotation types.

    The {annot_type} wildcard is either 'tsv_annotation' or 'refined_annotation'.
    """
    input:
        adata=f"{OUTDIR_PP}/{{sample}}/adata_{{sample}}_annotated.h5ad",
    output:
        results_dir=directory(f"{OUTDIR_PP}/{{sample}}/neighbourhood_analysis/{{annot_type}}"),
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
