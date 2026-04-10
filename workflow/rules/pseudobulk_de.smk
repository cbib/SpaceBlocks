rule pseudobulk_de:
    """
    Pseudobulk DE per cell type and region, run for both annotation types.

    The {annot_type} wildcard is either 'tsv_annotation' or 'refined_annotation',
    corresponding to obs columns 'cell_type_tsv' and 'cell_type_refined'.
    """
    input:
        adata=f"{OUTDIR_PP}/{{sample}}/adata_{{sample}}_annotated.h5ad",
    output:
        results_dir=directory(f"{OUTDIR_PP}/{{sample}}/pseudobulk_de/{{annot_type}}"),
    params:
        sample_id=lambda wc: wc.sample,
        annot_type=lambda wc: wc.annot_type,
    log:
        out=f"{LOGDIR}/pseudobulk_de/{{sample}}_{{annot_type}}.out",
        err=f"{LOGDIR}/pseudobulk_de/{{sample}}_{{annot_type}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_de/{{sample}}_{{annot_type}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("pseudobulk_de", "threads")
    resources:
        mem_mb=get_resource("pseudobulk_de", "mem_mb"),
        runtime=get_resource("pseudobulk_de", "runtime"),
    script:
        "../scripts/pseudobulk_de.py"
