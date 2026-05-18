rule ingest_ref:
    """Transfer cell-type labels from reference via sc.tl.ingest."""
    input:
        adata=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
        ingest_ref=config.get("ingest_ref", ""),
    output:
        adata_ingested=f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}_ingested.h5ad",
        plots_dir=directory(f"{SAMPLES_DIR}/{{sample}}/ingest"),
    params:
        sample_id=lambda wc: wc.sample,
        ref_label_key=config.get("ingest_ref_label_key", "cell_type"),
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
    log:
        out=f"{LOGDIR}/ingest_ref/{{sample}}.out",
        err=f"{LOGDIR}/ingest_ref/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/ingest_ref/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("ingest_ref", "threads")
    resources:
        mem_mb=get_resource("ingest_ref", "mem_mb"),
        runtime=get_resource("ingest_ref", "runtime"),
    script:
        "../scripts/ingest_ref.py"
