# workflow/rules/prepare_input_x5k.smk
# Xenium5k Headblock — zarr + QuPath GeoJSON → the standardized UNFILTERED CONTRACT
# h5ad (raw counts, embedded morphology image, region_annotation). No QC / normalise /
# cluster: that is the Coreblock's job (preprocess_umap). Analogous to prepare_input
# for Visium HD; writes the same single contract path every core consumer reads.
rule prepare_input_x5k:
    input:
        zarr=rules.convert_to_zarr.output.zarr,
        qupath_meta=rules.export_qupath_x5k.output.qupath_meta,   # geojson px→µm scale
    output:
        h5ad=config["contract"]["unfiltered_h5ad"],
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, input: str(input.zarr),
        geojson_dir=GEOJ_DIR,
        hires_pyramid_level=XENIUM5K.get("hires_pyramid_level", 3),
        pixel_size_um=XENIUM5K.get("pixel_size_um", 0.2125),
    log:
        out=f"{LOGDIR}/prepare_input_x5k/{{sample}}.out",
        err=f"{LOGDIR}/prepare_input_x5k/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/prepare_input_x5k/{{sample}}.tsv"
    conda:
        "../envs/xenium5k.yaml"
    threads:
        get_resource("prepare_input_x5k", "threads")
    resources:
        mem_mb=mem_mb_attempt("prepare_input_x5k"),
        runtime=get_resource("prepare_input_x5k", "runtime"),
    script:
        "../scripts/prepare_input_x5k.py"
