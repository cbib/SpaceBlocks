rule generate_qupath_vhd:
    """Copy the Space Ranger hires tissue image out for QuPath annotation."""
    input:
        sr_done=f"{OUTDIR_SR}/{{sample}}/.done",
    output:
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_tissue_hires_image.png",
    params:
        sample_id=lambda wc: wc.sample,
        sr_outdir=lambda wc, input: os.path.dirname(input.sr_done),
    log:
        out=f"{LOGDIR}/generate_qupath_vhd/{{sample}}.out",
        err=f"{LOGDIR}/generate_qupath_vhd/{{sample}}.err",
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("generate_qupath_vhd", "threads")
    resources:
        mem_mb=mem_mb_attempt("generate_qupath_vhd"),
        runtime=get_resource("generate_qupath_vhd", "runtime"),
    script:
        "../scripts/generate_qupath_vhd.py"
