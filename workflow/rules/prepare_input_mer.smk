rule prepare_input_mer:
    """
    MERSCOPE Headblock — Vizgen cell table + QuPath GeoJSON → the standardized UNFILTERED
    CONTRACT h5ad. Analogous to prepare_input_x5k; writes one AnnData.h5ad contract per
    sample. Depends on generate_qupath_mer for the coordinate metadata and grey background
    (so the mosaic is read once and the contract shares the QuPath grid).
    """
    input:
        qupath_meta=rules.generate_qupath_mer.output.qupath_meta,   # px<->µm mapping + p0
        background=rules.generate_qupath_mer.output.background,      # grey contract image
    output:
        h5ad=config["contract"]["unfiltered_h5ad"],
    params:
        sample_id=lambda wc: wc.sample,
        merscope_dir=lambda wc: merscope_dir_for(wc.sample),
        geojson_dir=GEOJ_DIR,
    log:
        out=f"{LOGDIR}/prepare_input_mer/{{sample}}.out",
        err=f"{LOGDIR}/prepare_input_mer/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/prepare_input_mer/{{sample}}.tsv"
    conda:
        "../envs/xenium5k.yaml"
    threads:
        get_resource("prepare_input_mer", "threads")
    resources:
        mem_mb=mem_mb_attempt("prepare_input_mer"),
        runtime=get_resource("prepare_input_mer", "runtime"),
    script:
        "../scripts/prepare_input_mer.py"
