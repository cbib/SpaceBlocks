rule process_geojson:
    """Load a sample's OPTIONAL region GeoJSON ({sample}_tissue_hires_image.geojson
    or {sample}_morphology.geojson) and write it into the sample dir. Only invoked
    for samples in SAMPLES_WITH_GEOJSON (see common.smk / _find_geojson)."""
    input:
        geojson=lambda wc: _find_geojson(wc.sample),
    output:
        geojson=f"{SAMPLES_DIR}/{{sample}}/{{sample}}_annotations.geojson",
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