rule generate_qupath_ate:
    """
    Atera Headblock — composite a morphology RGB TIFF (+ scale-factor JSON) from the
    zarr for QuPath region annotation. This is the DEFAULT annotation image and is
    always produced; the optional registered H&E is handled by generate_qupath_he_ate.
    Analogous to generate_qupath_x5k (Xenium 5K).
    """
    input:
        done=rules.convert_zarr_ate.output.done,
    output:
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology.tiff",
        qupath_meta=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology_scalefactors.json",
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, input: str(input.done)[:-len(".done")],
        qupath_pyramid_level=ATERA.get("qupath_pyramid_level", 3),
        pixel_size_um=ATERA.get("pixel_size_um", 0.2125),
    log:
        out=f"{LOGDIR}/generate_qupath_ate/{{sample}}.out",
        err=f"{LOGDIR}/generate_qupath_ate/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/generate_qupath_ate/{{sample}}.tsv"
    conda:
        "../envs/atera.yaml"
    threads:
        get_resource("generate_qupath_ate", "threads")
    resources:
        mem_mb=mem_mb_attempt("generate_qupath_ate"),
        runtime=get_resource("generate_qupath_ate", "runtime"),
    script:
        "../scripts/generate_qupath_ate.py"
