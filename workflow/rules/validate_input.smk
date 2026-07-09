# workflow/rules/validate_input.smk
# ─────────────────────────────────────────────────────────────────────────────
# Input-contract validation for the modular CORE. Run once per sample, BEFORE
# qc_sweep and run_preprocessing. It asserts the unfiltered, segmented h5ad satisfies
# the core contract and records whether an image is embedded (so spatial plots
# can later choose sc.pl.spatial-with-image vs a coordinate scatter). It does NOT
# modify the data — file preparation is a HEAD concern.
#
# If validation fails, the script raises and no report is written, so Snakemake
# aborts the DAG before any analysis runs. On success it writes a JSON report
# whose `image_present` / `image_mode` fields downstream rules can read.
#
# Expected config (example — wire these into your top-level config.yaml):
#
#   outdir: "results"
#   samples: [S1, S2, S3]                 # however you enumerate samples
#   contract:
#     unfiltered_h5ad: "results/head/{sample}/unfiltered.h5ad"  # {sample} pattern
#     sample_key: "sample"                # obs column identifying the sample
#     spatial_key: "spatial"              # obsm key for coords / uns key for image
#     require_region: false               # error if region_annotation absent
#     require_raw_counts: true            # error if X is not integer-valued
#     mito_prefix: ["MT-", "mt-"]         # gene-name prefixes for mito detection
# ─────────────────────────────────────────────────────────────────────────────

_CONTRACT = config.get("contract", {})
_OUT = config.get("outdir", "results")


rule validate_input:
    """Assert the unfiltered h5ad satisfies the core contract; record image presence."""
    input:
        h5ad=lambda wc: _CONTRACT["unfiltered_h5ad"].format(sample=wc.sample),
    output:
        report=f"{_OUT}/{{sample}}/validation/input_validation.json",
    params:
        sample_id=lambda wc: wc.sample,
        sample_key=_CONTRACT.get("sample_key", "sample"),
        spatial_key=_CONTRACT.get("spatial_key", "spatial"),
        require_region=_CONTRACT.get("require_region", False),
        require_raw_counts=_CONTRACT.get("require_raw_counts", True),
        mito_prefix=_CONTRACT.get("mito_prefix", ["MT-", "mt-"]),
    log:
        out=f"{_OUT}/logs/{{sample}}/validate_input.log",
        err=f"{_OUT}/logs/{{sample}}/validate_input.err",
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("validate_input", "threads")
    resources:
        mem_mb=mem_mb_attempt("validate_input"),
        runtime=get_resource("validate_input", "runtime"),
    script:
        "../scripts/validate_input.py"
