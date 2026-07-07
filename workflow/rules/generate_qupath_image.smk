# workflow/rules/generate_qupath_image.smk
# ─────────────────────────────────────────────────────────────────────────────
# Visium HD HEAD (step 2). Copies the Space Ranger hires tissue image out of the
# (version-varying) SR tree to a stable path, so it can be opened and annotated
# in QuPath. The exported "{sample}_tissue_hires_image.geojson" (placed in
# GEOJ_DIR) is later consumed by prepare_input.
#
# Kept SEPARATE from prepare_input on purpose: it is a cheap file copy (low
# resources) and is meant to be grouped into a single DAG job together with
# spaceranger — the geojson annotation must exist BEFORE prepare_input runs.
# (DAG grouping will be wired in a later revision.)
# Globals used (defined in the Snakefile): OUTDIR_SR, SAMPLES_DIR, LOGDIR, get_resource.
# ─────────────────────────────────────────────────────────────────────────────

rule generate_qupath_image:
    """Copy the Space Ranger hires tissue image out for QuPath annotation."""
    input:
        sr_done=f"{OUTDIR_SR}/{{sample}}/.done",
    output:
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_tissue_hires_image.png",
    params:
        sample_id=lambda wc: wc.sample,
        sr_outdir=lambda wc, input: os.path.dirname(input.sr_done),
    log:
        out=f"{LOGDIR}/generate_qupath_image/{{sample}}.out",
        err=f"{LOGDIR}/generate_qupath_image/{{sample}}.err",
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("generate_qupath_image", "threads")
    resources:
        mem_mb=get_resource("generate_qupath_image", "mem_mb"),
        runtime=get_resource("generate_qupath_image", "runtime"),
    script:
        "../scripts/generate_qupath_image.py"
