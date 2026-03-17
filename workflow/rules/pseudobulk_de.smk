rule pseudobulk_de:
    """
    Pseudobulk differential expression per cell type and region.

    Requires a user-provided cluster annotation TSV mapping Leiden clusters
    to cell-type labels.  Aggregates raw counts per cell_type × sample
    (and × region_annotation if available), then runs pyDESeq2.

    This rule only executes if the annotation file exists and is non-empty;
    the input function in the Snakefile handles this gating.
    """
    input:
        adata=f"{OUTDIR_PP}/{{sample}}/adata_{{sample}}.h5ad",
        cluster_annotations=config.get("cluster_annotations", ""),
    output:
        results_dir=directory(f"{OUTDIR_PP}/{{sample}}/pseudobulk_de"),
    params:
        sample_id=lambda wc: wc.sample,
    log:
        out=f"{LOGDIR}/pseudobulk_de/{{sample}}.out",
        err=f"{LOGDIR}/pseudobulk_de/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/pseudobulk_de/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("pseudobulk_de", "threads")
    resources:
        mem_mb=get_resource("pseudobulk_de", "mem_mb"),
        runtime=get_resource("pseudobulk_de", "runtime"),
    script:
        "../scripts/pseudobulk_de.py"
