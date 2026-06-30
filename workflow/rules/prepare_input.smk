# workflow/rules/prepare_input.smk
# ─────────────────────────────────────────────────────────────────────────────
# Visium HD HEAD (step 3). Builds the UNFILTERED contract h5ad the common core
# consumes (raw counts + obsm["spatial"] + uns image + region_annotation) from
# Space Ranger outputs and the QuPath geojson annotations.
#
# Runs AFTER generate_qupath_image + manual QuPath annotation (geojson in
# GEOJ_DIR). The output path is config["contract"]["unfiltered_h5ad"] (set in the
# Snakefile), so prepare_input writes exactly what validate_input / qc_sweep /
# preprocess_umap read. Globals used (Snakefile): OUTDIR_SR, GEOJ_DIR, LOGDIR,
# get_resource.
# ─────────────────────────────────────────────────────────────────────────────

rule prepare_input:
    """Build the UNFILTERED contract h5ad from Space Ranger + QuPath annotations."""
    input:
        sr_done=f"{OUTDIR_SR}/{{sample}}/.done",
    output:
        h5ad=config["contract"]["unfiltered_h5ad"],
    params:
        sample_id=lambda wc: wc.sample,
        sr_outdir=f"{OUTDIR_SR}/{{sample}}",
        geojson_path=GEOJ_DIR,
    log:
        out=f"{LOGDIR}/prepare_input/{{sample}}.out",
        err=f"{LOGDIR}/prepare_input/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/prepare_input/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("prepare_input", "threads")
    resources:
        mem_mb=get_resource("prepare_input", "mem_mb"),
        runtime=get_resource("prepare_input", "runtime"),
    script:
        "../scripts/prepare_input.py"
