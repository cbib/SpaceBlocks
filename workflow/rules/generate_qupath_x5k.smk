rule generate_qupath_x5k:
    """
    Xenium5k Headblock — composite a morphology RGB TIFF (+ scale-factor JSON) from the
    zarr for QuPath region annotation. Analogous to generate_qupath_vhd (Visium HD).
    """
    input:
        done=rules.convert_zarr_x5k.output.done,
    output:
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology.tiff",
        qupath_meta=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology_scalefactors.json",
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, input: str(input.done)[:-len(".done")],
        qupath_pyramid_level=XENIUM5K.get("qupath_pyramid_level", 3),
    log:
        out=f"{LOGDIR}/generate_qupath_x5k/{{sample}}.out",
        err=f"{LOGDIR}/generate_qupath_x5k/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/generate_qupath_x5k/{{sample}}.tsv"
    conda:
        "../envs/xenium5k.yaml"
    threads:
        get_resource("generate_qupath_x5k", "threads")
    resources:
        mem_mb=mem_mb_attempt("generate_qupath_x5k"),
        runtime=get_resource("generate_qupath_x5k", "runtime"),
    script:
        "../scripts/generate_qupath_x5k.py"
