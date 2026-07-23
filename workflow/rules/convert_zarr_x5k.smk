# workflow/rules/convert_zarr_x5k.smk
# Xenium5k Headblock — heavy I/O: convert one Xenium output bundle to a SpatialData
# zarr store (the canonical object holding morphology image, masks, boundaries, and
# the raw count table). Analogous to spaceranger_count_vhd for Visium HD.
rule convert_zarr_x5k:
    input:
        xenium_dir=lambda wc: xenium_dir_for(wc.sample),
    output:
        # A .done marker (touched by the script only after sdata.write() succeeds) is the
        # tracked output — the zarr store itself is a directory the script manages (it
        # cleans a partial one on re-run), so completion is precise and version-agnostic.
        done=f"{XENIUM_ZARR_DIR}/{{sample}}/{{sample}}.zarr.done",
    params:
        sample_id=lambda wc: wc.sample,
        zarr_path=lambda wc, output: str(output.done)[:-len(".done")],
    log:
        out=f"{LOGDIR}/convert_zarr_x5k/{{sample}}.out",
        err=f"{LOGDIR}/convert_zarr_x5k/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/convert_zarr_x5k/{{sample}}.tsv"
    conda:
        "../envs/xenium5k.yaml"
    threads:
        get_resource("convert_zarr_x5k", "threads")
    resources:
        mem_mb=mem_mb_attempt("convert_zarr_x5k"),
        runtime=get_resource("convert_zarr_x5k", "runtime"),
    script:
        "../scripts/convert_zarr_x5k.py"
