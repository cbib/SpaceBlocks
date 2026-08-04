rule convert_zarr_ate:
    """
    Atera Headblock — heavy I/O: convert one Atera output bundle to a SpatialData
    zarr store. Analogous to convert_zarr_x5k (Xenium 5K) and spaceranger_count_vhd
    (Visium HD).
    """
    input:
        atera_dir=lambda wc: atera_dir_for(wc.sample),
    output:
        done=f"{ATERA_ZARR_DIR}/{{sample}}/{{sample}}.zarr.done",
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, output: str(output.done)[:-len(".done")],
    log:
        out=f"{LOGDIR}/convert_zarr_ate/{{sample}}.out",
        err=f"{LOGDIR}/convert_zarr_ate/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/convert_zarr_ate/{{sample}}.tsv"
    conda:
        "../envs/atera.yaml"
    threads:
        get_resource("convert_zarr_ate", "threads")
    resources:
        mem_mb=mem_mb_attempt("convert_zarr_ate"),
        runtime=get_resource("convert_zarr_ate", "runtime"),
    script:
        "../scripts/convert_zarr_ate.py"
