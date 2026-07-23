# workflow/rules/prepare_input_vhd.smk
# ─────────────────────────────────────────────────────────────────────────────
# Visium HD HEAD (step 3). Builds the UNFILTERED contract h5ad the common core
# consumes (raw counts + obsm["spatial"] + uns image + region_annotation) from
# Space Ranger outputs and the QuPath geojson annotations.
#
# Runs AFTER generate_qupath_vhd + manual QuPath annotation (geojson in
# GEOJ_DIR). The output path is config["contract"]["unfiltered_h5ad"] (set in the
# Snakefile), so prepare_input_vhd writes exactly what validate_input / qc_sweep /
# preprocess_umap read. Globals used (Snakefile): OUTDIR_SR, GEOJ_DIR, LOGDIR,
# get_resource.
# ─────────────────────────────────────────────────────────────────────────────

rule prepare_input_vhd:
    """Build the UNFILTERED contract h5ad from Space Ranger + QuPath annotations."""
    input:
        sr_done=f"{OUTDIR_SR}/{{sample}}/.done",
    output:
        h5ad=config["contract"]["unfiltered_h5ad"],
    params:
        sample_id=lambda wc: wc.sample,
        sr_outdir=lambda wc, input: os.path.dirname(input.sr_done),
        geojson_path=GEOJ_DIR,
    log:
        out=f"{LOGDIR}/prepare_input_vhd/{{sample}}.out",
        err=f"{LOGDIR}/prepare_input_vhd/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/prepare_input_vhd/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("prepare_input_vhd", "threads")
    resources:
        mem_mb=mem_mb_attempt("prepare_input_vhd"),
        runtime=get_resource("prepare_input_vhd", "runtime"),
    script:
        "../scripts/prepare_input_vhd.py"
