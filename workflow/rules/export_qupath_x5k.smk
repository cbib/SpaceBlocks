# workflow/rules/export_qupath_x5k.smk
# Xenium5k Headblock — composite a morphology RGB TIFF (+ scale-factor JSON) from the
# zarr for QuPath region annotation. Analogous to generate_qupath_image (Visium HD):
# run it, annotate in QuPath, export {sample}_morphology.geojson into geojson_path,
# then continue to prepare_input_x5k.
rule export_qupath_x5k:
    input:
        zarr=rules.convert_to_zarr.output.zarr,
    output:
        qupath_image=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology.tiff",
        qupath_meta=f"{SAMPLES_DIR}/{{sample}}/QuPath_image/{{sample}}_morphology_scalefactors.json",
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, input: str(input.zarr),
        qupath_pyramid_level=XENIUM5K.get("qupath_pyramid_level", 3),
    log:
        out=f"{LOGDIR}/export_qupath_x5k/{{sample}}.out",
        err=f"{LOGDIR}/export_qupath_x5k/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/export_qupath_x5k/{{sample}}.tsv"
    conda:
        "../envs/xenium5k.yaml"
    threads:
        get_resource("export_qupath_x5k", "threads")
    resources:
        mem_mb=mem_mb_attempt("export_qupath_x5k"),
        runtime=get_resource("export_qupath_x5k", "runtime"),
    script:
        "../scripts/export_qupath_x5k.py"
