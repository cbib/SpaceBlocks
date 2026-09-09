rule process_geojson:
    """Validate an optional region GeoJSON without modifying the source file."""
    input:
        geojson=lambda wc: _find_geojson(wc.sample),
    output:
        report=f"{SAMPLES_DIR}/{{sample}}/validation/geojson_validation.json",
    log:
        out=f"{LOGDIR}/process_geojson/{{sample}}.out",
        err=f"{LOGDIR}/process_geojson/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/process_geojson/{{sample}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads: get_resource("process_geojson", "threads")
    resources:
        mem_mb=mem_mb_attempt("process_geojson"),
        runtime=get_resource("process_geojson", "runtime"),
    params:
        sample_id=lambda wc: wc.sample,
    script:
        "../scripts/process_geojson.py"
