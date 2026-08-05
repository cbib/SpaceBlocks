rule generate_qupath_mer:
    """
    MERSCOPE Headblock — composite a morphology RGB TIFF (+ pixel<->micron JSON and a grey
    background reused by the contract) from the Vizgen mosaic OME-TIFFs for QuPath region
    annotation. Analogous to generate_qupath_x5k, but reads the mosaic directly (no zarr).
    """
    input:
        # anchor on the always-present transform CSV; the script discovers the mosaics
        transform=lambda wc: os.path.join(
            merscope_dir_for(wc.sample), "images", "micron_to_mosaic_pixel_transform.csv"),
    output:
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology.tiff",
        qupath_meta=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology_scalefactors.json",
        background=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_background.npy",
    params:
        sample_id=lambda wc: wc.sample,
        merscope_dir=lambda wc: merscope_dir_for(wc.sample),
        z_index=MERSCOPE.get("z_index", 3),
        hires_pixel_size_um=MERSCOPE.get("hires_pixel_size_um", 1.0),
        channels=MERSCOPE.get("channels", DEFAULT_MER_CHANNELS),
    log:
        out=f"{LOGDIR}/generate_qupath_mer/{{sample}}.out",
        err=f"{LOGDIR}/generate_qupath_mer/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/generate_qupath_mer/{{sample}}.tsv"
    conda:
        "../envs/xenium5k.yaml"
    threads:
        get_resource("generate_qupath_mer", "threads")
    resources:
        mem_mb=mem_mb_attempt("generate_qupath_mer"),
        runtime=get_resource("generate_qupath_mer", "runtime"),
    script:
        "../scripts/generate_qupath_mer.py"
