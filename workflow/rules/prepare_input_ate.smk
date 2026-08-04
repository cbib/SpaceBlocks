rule prepare_input_ate:
    """
    Atera Headblock — zarr + QuPath GeoJSON → the standardized UNFILTERED CONTRACT
    h5ad. Analogous to prepare_input_x5k (Xenium 5K) and prepare_input_vhd (Visium HD);
    writes one AnnData .h5ad contract file per sample.

    Inputs are assembled by _ate_prepare_inputs so that the H&E scale-factor JSON is
    required only when the optional H&E rule is active.
    """
    input:
        unpack(_ate_prepare_inputs),
    output:
        h5ad=config["contract"]["unfiltered_h5ad"],
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, input: str(input.done)[:-len(".done")],
        geojson_dir=GEOJ_DIR,
        hires_pyramid_level=ATERA.get("hires_pyramid_level", 3),
        pixel_size_um=ATERA.get("pixel_size_um", 0.2125),
    log:
        out=f"{LOGDIR}/prepare_input_ate/{{sample}}.out",
        err=f"{LOGDIR}/prepare_input_ate/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/prepare_input_ate/{{sample}}.tsv"
    conda:
        "../envs/atera.yaml"
    threads:
        get_resource("prepare_input_ate", "threads")
    resources:
        mem_mb=mem_mb_attempt("prepare_input_ate"),
        runtime=get_resource("prepare_input_ate", "runtime"),
    script:
        "../scripts/prepare_input_ate.py"
