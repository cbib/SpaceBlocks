_CONTRACT = config.get("contract", {})
_OUT = config.get("outdir", "results")
_QC = config.get("qc_sweep", {})


rule qc_sweep:
    """Per-sample QC diagnostic: violins (per region), spatial threshold-landing
    grid, continuous QC-feature maps, joint scatter, and threshold tables."""
    input:
        h5ad=lambda wc: _CONTRACT["unfiltered_h5ad"].format(sample=wc.sample),
        validation=f"{_OUT}/{{sample}}/validation/input_validation.json",
        ingest_ref=lambda wc: (config.get("ingest_ref", "")
                               if _QC.get("ingest_enabled") and config.get("ingest_ref")
                               else []),
    output:
        violins_png=f"{_OUT}/{{sample}}/qc_sweep/qc_violins.png",
        low_genes_png=f"{_OUT}/{{sample}}/qc_sweep/qc_remove_low_genes.png",
        low_counts_png=f"{_OUT}/{{sample}}/qc_sweep/qc_remove_low_counts.png",
        high_counts_png=f"{_OUT}/{{sample}}/qc_sweep/qc_remove_high_counts.png",
        high_mito_png=f"{_OUT}/{{sample}}/qc_sweep/qc_remove_high_mito.png",
        features_png=f"{_OUT}/{{sample}}/qc_sweep/qc_features_spatial.png",
        joint_png=f"{_OUT}/{{sample}}/qc_sweep/qc_joint_scatter.png",
        summary_tsv=f"{_OUT}/{{sample}}/qc_sweep/qc_thresholds_summary.tsv",
        ingest_tsv=f"{_OUT}/{{sample}}/qc_sweep/qc_ingest_removed.tsv",
    params:
        sample_id=lambda wc: wc.sample,
        thresholds=lambda wc: _thresholds_for(wc.sample),
        mito_prefix=_CONTRACT.get("mito_prefix", ["MT-", "mt-"]),
        dpi=_QC.get("dpi", config.get("analysis", {}).get("plot_dpi", 300)),
        ingest_enabled=_QC.get("ingest_enabled", False),
        ref_label_key=config.get("ingest_ref_label_key", "cell_type"),
        region_colors=config.get("analysis", {}).get("region_colors", {}),
        region_levels=config.get("analysis", {}).get("region_levels", []),
    log:
        out=f"{_OUT}/logs/{{sample}}/qc_sweep.log",
        err=f"{_OUT}/logs/{{sample}}/qc_sweep.err",
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("qc_sweep", "threads")
    resources:
        mem_mb=mem_mb_attempt("qc_sweep"),
        runtime=get_resource("qc_sweep", "runtime"),
    script:
        "../scripts/qc_sweep.py"
