rule generate_annotation_template:
    """
    Generate an empty cluster annotation template TSV.

    Rows: cluster numbers 0–49.
    Columns: {sample_id}_{resolution} for every sample × resolution
    combination in the scan range.

    The user fills in cell-type labels, then sets ``cluster_annotations``
    in config.yaml to this file to trigger pseudobulk_de and
    neighbourhood_analysis.
    """
    output:
        template=f"{OUTDIR_PP}/cluster_annotations_template.tsv",
    params:
        sample_ids=SAMPLE_IDS,
        resolutions=RESOLUTIONS,
        max_clusters=50,
    log:
        out=f"{LOGDIR}/generate_annotation_template.out",
        err=f"{LOGDIR}/generate_annotation_template.err",
    conda:
        "../envs/visiumhd.yaml"
    localrule: True
    script:
        "../scripts/generate_annotation_template.py"
