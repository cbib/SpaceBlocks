SUBCOMPARTMENTS = list(config.get("subcompartments", {}).keys())

rule subcluster:
    """
    Subset and subcluster a cell compartment from the concatenated dataset.

    Produces two branches (Harmony / NoHarmony), each with leiden clustering
    at multiple resolutions, silhouette evaluation, QC plots, and UMAPs.
    """
    input:
        adata=f"{OUTDIR_PP}/integrated_samples/concatenated.h5ad",
    output:
        sub_dir=directory(f"{OUTDIR_PP}/Subcompartments/{{subcompartment}}"),
    wildcard_constraints:
        subcompartment="|".join(SUBCOMPARTMENTS) if SUBCOMPARTMENTS else "NONE",
    params:
        subcompartment=lambda wc: wc.subcompartment,
        strings=lambda wc: config["subcompartments"][wc.subcompartment]["strings"],
        annot_col=config.get("subcompartment_annot_col", "cell_type_tsv"),
        resolution_min=lambda wc: config["subcompartments"][wc.subcompartment].get("resolution_min", 0.2),
        resolution_max=lambda wc: config["subcompartments"][wc.subcompartment].get("resolution_max", 1.0),
        resolution_step=lambda wc: config["subcompartments"][wc.subcompartment].get("resolution_step", 0.2),
        n_neighbors=ANALYSIS.get("n_neighbors", 10),
        n_pcs=ANALYSIS.get("n_pcs", 30),
        de_n_genes=ANALYSIS.get("de_n_genes", 10),
        random_seed=RANDOM_SEED,
        annotation_colors=config.get("annotation_colors", {}),
        region_colors=ANALYSIS.get("region_colors", {}),
        niche_column=GENE_EXPLORATION.get("niche_column", ""),
    log:
        out=f"{LOGDIR}/subcluster/{{subcompartment}}.out",
        err=f"{LOGDIR}/subcluster/{{subcompartment}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/subcluster/{{subcompartment}}.tsv"
    conda:
        "../envs/visiumhd.yaml"
    threads:
        get_resource("subcluster", "threads")
    resources:
        mem_mb=get_resource("subcluster", "mem_mb"),
        runtime=get_resource("subcluster", "runtime"),
    script:
        "../scripts/subcluster.py"
