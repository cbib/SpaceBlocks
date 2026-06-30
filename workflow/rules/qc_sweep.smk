# workflow/rules/qc_sweep.smk
# ─────────────────────────────────────────────────────────────────────────────
# OPTIONAL pre-filtering QC diagnostic (modular core, module 1). Per sample, it
# reads the VALIDATED unfiltered h5ad and shows where candidate QC cutoffs land
# on the tissue, so the final analysis.* thresholds can be chosen deliberately.
# It does NOT filter, NOT cluster and does NOT write an h5ad. Request its outputs
# explicitly — it is not pulled in by run_upstream.
#
# Thresholds are LISTS per feature and may be SHARED across samples (a `default`
# block) or UNIQUE per sample (a `per_sample` override block, merged key-by-key).
# An ingest reference (cell-type labels) can optionally be overlaid for a richer
# "what am I removing" diagnosis.
#
# Expected config (in addition to the `contract` / `outdir` / `samples` blocks
# used by validate_input.smk):
#
#   qc_sweep:
#     ingest_enabled: false
#     ingest_pattern: "results/head/{sample}/ingested.h5ad"   # {sample} pattern
#     ingest_key: "cell_type_ingest"
#     cell_id_key: "cell_id"
#     dpi: 300
#     thresholds:
#       default:
#         min_genes:  [100, 200]
#         min_counts: [10, 50]
#         max_counts: [3000, 5000]
#         max_pct_mt: [10, 15, 20]
#       per_sample:                      # optional, overrides default key-by-key
#         SAMPLE_A:
#           max_counts: [4000]
# ─────────────────────────────────────────────────────────────────────────────

_CONTRACT = config.get("contract", {})
_OUT = config.get("outdir", "results")
_QC = config.get("qc_sweep", {})
_TH = _QC.get("thresholds", {})


def _thresholds_for(sample):
    """Shared `default` thresholds overridden key-by-key by any `per_sample` set."""
    merged = dict(_TH.get("default", {}))
    merged.update(_TH.get("per_sample", {}).get(sample, {}))
    return merged


rule qc_sweep:
    """Per-sample QC diagnostic: violins (per region), spatial threshold-landing
    grid, continuous QC-feature maps, joint scatter, and threshold tables."""
    input:
        h5ad=lambda wc: _CONTRACT["unfiltered_h5ad"].format(sample=wc.sample),
        validation=f"{_OUT}/{{sample}}/validation/input_validation.json",
        ingest=lambda wc: (_QC["ingest_pattern"].format(sample=wc.sample)
                           if _QC.get("ingest_enabled") else []),
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
        ingest_key=_QC.get("ingest_key", "cell_type_ingest"),
        cell_id_key=_QC.get("cell_id_key", "cell_id"),
    log:
        out=f"{_OUT}/logs/{{sample}}/qc_sweep.log",
        err=f"{_OUT}/logs/{{sample}}/qc_sweep.err",
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("qc_sweep", "threads")
    resources:
        mem_mb=get_resource("qc_sweep", "mem_mb"),
        runtime=get_resource("qc_sweep", "runtime"),
    script:
        "../scripts/qc_sweep.py"
