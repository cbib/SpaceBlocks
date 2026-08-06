rule generate_qupath_he_ate:
    """
    Atera Headblock (OPTIONAL) — downsample the registered H&E whole-slide image into a
    QuPath-annotatable TIFF, compose the H&E-pixel → micron transform from the 10x
    alignment matrix, and resample the H&E onto the morphology pixel grid to serve as
    the contract background image.

    Only included when atera.he_image and atera.he_alignment are configured. The H&E
    ships as a separate download from the outs/ bundle, so its paths are given
    explicitly in the config (absolute {sample} patterns) rather than discovered.

    Deliberately standalone rather than folded into generate_qupath_ate: the H&E is
    optional (a conditional output block is not expressible cleanly in Snakemake), and
    keeping it separate makes it liftable into the Xenium head, which supports the same
    aligned-image workflow.
    """
    input:
        done=rules.convert_zarr_ate.output.done,
        he_image=lambda wc: he_file_for(wc.sample, "he_image"),
        he_alignment=lambda wc: he_file_for(wc.sample, "he_alignment"),
    output:
        qupath_meta=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_he_scalefactors.json",
        he_background=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_he_background.tiff",
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, input: str(input.done)[:-len(".done")],
        # Optional: not an input: — absence must not block the rule, it only disables QA.
        he_keypoints=lambda wc: he_file_for(wc.sample, "he_keypoints", required=False),
        he_pyramid_level=ATERA.get("he_pyramid_level", 4),
        # Must match prepare_input_ate's level: the background is warped onto that grid.
        hires_pyramid_level=ATERA.get("hires_pyramid_level", 3),
        pixel_size_um=ATERA.get("pixel_size_um", 0.2125),
        residual_warn_px=ATERA.get("he_residual_warn_px", 50),
    log:
        out=f"{LOGDIR}/generate_qupath_he_ate/{{sample}}.out",
        err=f"{LOGDIR}/generate_qupath_he_ate/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/generate_qupath_he_ate/{{sample}}.tsv"
    conda:
        "../envs/atera.yaml"
    threads:
        get_resource("generate_qupath_he_ate", "threads")
    resources:
        mem_mb=mem_mb_attempt("generate_qupath_he_ate"),
        runtime=get_resource("generate_qupath_he_ate", "runtime"),
    script:
        "../scripts/generate_qupath_he_ate.py"
