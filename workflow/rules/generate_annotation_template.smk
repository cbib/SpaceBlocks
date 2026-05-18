rule generate_annotation_template:
    """Generate empty cluster annotation template TSV."""
    output:
        template=f"{OUTDIR_PP}/cluster_annotations_template.tsv",
    params:
        sample_ids=SAMPLE_IDS,
        resolutions=RESOLUTIONS,
        max_clusters=50,
    log:
        out=f"{LOGDIR}/generate_annotation_template.out",
        err=f"{LOGDIR}/generate_annotation_template.err",
    localrule: True
    run:
        import pandas as pd
        columns = [f"{s}_{r}" for s in params.sample_ids for r in params.resolutions]
        df = pd.DataFrame(index=range(params.max_clusters), columns=columns, data="")
        df.index.name = "cluster"
        df.to_csv(str(output.template), sep="\t")
