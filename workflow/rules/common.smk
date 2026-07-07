# workflow/rules/common.smk
# ─────────────────────────────────────────────────────────────────────────────
# Shared helper functions for the whole workflow. Collected here (rather than mixed
# in with rules) per the Snakemake modularization guideline: rule files contain
# only rules, this file contains only functions. Included FIRST from the Snakefile,
# so every function is visible to the globals block and to all rule modules. The
# functions reference workflow globals (config, SAMPLE_IDS, SAMPLES_DIR, …) that the
# Snakefile defines after this include; that is fine because the functions are only
# *called* later (at global-evaluation or DAG-build time), never at include time.
# ─────────────────────────────────────────────────────────────────────────────


# ── Resources ────────────────────────────────────────────────────────────────
def get_resource(rule_name, key):
    """Per-rule resource lookup with a `default` fallback (config['resources'])."""
    section = config.get("resources", {})
    rule_res = section.get(rule_name, section.get("default", {}))
    if key not in rule_res:
        sys.exit(f"[config error] Resource '{key}' missing for rule '{rule_name}'.")
    return rule_res[key]


# ── Small derivations called from the Snakefile globals block ────────────────
def _resolution_range(rmin, rmax, step):
    """Inclusive Leiden resolution ladder from (min, max, step)."""
    n = round((rmax - rmin) / step)
    return [round(rmin + i * step, 10) for i in range(n + 1)]


def _cluster_annotations_available():
    """True when config['cluster_annotations'] points at a non-empty TSV."""
    path = config.get("cluster_annotations", "")
    if not path or not os.path.isfile(path):
        return False
    try:
        df = pd.read_csv(path, sep="\t", index_col=0)
        return df.shape[0] > 0 and df.shape[1] > 0 and not df.isna().all().all()
    except Exception:
        return False


def _extra_annot_color_triples():
    """Flatten sample_colors for the extra-annotation columns into three parallel
    lists (column, value, hex) for safe Python→R hand-off in pseudobulk_de."""
    cols, vals, colors = [], [], []
    sc_ = SAMPLE_COLORS if isinstance(SAMPLE_COLORS, dict) else {}
    for c in EXTRA_ANNOT_COLUMNS:
        for v, hexc in (sc_.get(c, {}) or {}).items():
            cols.append(c); vals.append(str(v)); colors.append(str(hexc))
    return cols, vals, colors


def core_sample_meta(sample):
    """Experimental-design metadata for `sample` (design columns only, excludes the
    `sample` key). Returns {} when no core sheet is configured."""
    if CORE_SAMPLES is None or sample not in CORE_SAMPLES.index:
        return {}
    row = CORE_SAMPLES.loc[sample].to_dict()
    row.pop("sample", None)
    return row


# ── HEAD (Visium HD) input helpers ───────────────────────────────────────────
def fastq_dirs_comma_separated(wc):
    """Comma-separated fastq dirs for spaceranger (base + any re-sequencing runs)."""
    s = SAMPLES[wc.sample]
    dirs = [f"{s['fastq_dir']}/fastqs/"]
    for col in ("resequenced_dir", "new_resequenced_dir", "third_resequencing_dir"):
        val = s.get(col, "")
        if val not in ("", "nan"):
            dirs.append(f"{val}/fastqs/")
    return ",".join(dirs)


# ── CORE per-rule input functions ────────────────────────────────────────────
def _thresholds_for(sample):
    """Global candidate threshold lists for qc_sweep (flat: feature -> list of
    cut-offs). One set is applied to every sample; per-sample QC tuning lives in
    preprocess_umap."""
    return dict(config.get("qc_sweep", {}).get("thresholds", {}))


def _preprocess_inputs(wc):
    """preprocess_umap inputs: the validated contract h5ad (gated on the validation
    report), plus OPTIONAL precomputed metadata and per-sample QC-threshold TSV,
    each declared only when present so Snakemake tracks edits without hard-requiring
    them."""
    inputs = {
        "h5ad": config["contract"]["unfiltered_h5ad"].format(sample=wc.sample),
        # gate: core runs only after the contract is validated
        "validation": f"{config['outdir']}/{wc.sample}/validation/input_validation.json",
    }
    if USE_PRECOMPUTED:
        precomp_dir = config.get("precomputed_metadata_dir", "")
        if precomp_dir:
            meta_file = os.path.join(precomp_dir, f"metadata_{wc.sample}.tsv")
            if os.path.isfile(meta_file):
                inputs["precomputed_meta"] = meta_file
    th_tsv = config.get("preprocess_thresholds", "") or ""
    if th_tsv and os.path.isfile(th_tsv):
        inputs["thresholds_tsv"] = th_tsv
    return inputs


def _spatial_niches_inputs(wc):
    """Inputs for spatial_niches. When use_precomputed is set with an external
    niche_dir, the per-sample niche_{sample}.tsv are added as MANDATORY inputs so
    Snakemake aborts before the job if any is absent (complementary to the in-script
    safeguard, which also catches empty/partial TSVs). When niche_dir is empty the
    rule reloads from its own output dir, so nothing is added here (no circular
    dependency)."""
    inputs = {
        "adatas": expand(f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad",
                         sample=SAMPLE_IDS),
    }
    sn = config.get("spatial_niches", {})
    if sn.get("use_precomputed", False):
        d = sn.get("niche_dir", "")
        if d:
            inputs["precomputed_niches"] = [
                os.path.join(d, f"niche_{s}.tsv") for s in SAMPLE_IDS
            ]
    return inputs


def _annotate_input_adata(wc):
    """annotate_cells reads the ingested adata when an ingest reference is
    configured, else the plain preprocessed adata."""
    if config.get("ingest_ref", ""):
        return f"{SAMPLES_DIR}/{wc.sample}/adata_{wc.sample}_ingested.h5ad"
    return f"{SAMPLES_DIR}/{wc.sample}/adata_{wc.sample}.h5ad"


def _annotate_niche_input(wc):
    """Spatial-niche TSV for this sample (barcode → spatial_niche), injected as the
    spatial_niche obs column. Always the spatial_niches rule output (the rule handles
    precomputed reload internally). Returns [] when niche identification is disabled,
    so annotate_cells does not depend on it."""
    sn = config.get("spatial_niches", {})
    if not sn.get("enabled", False):
        return []
    return f"{OUTDIR_PP}/spatial_niches/tsv/niche_{wc.sample}.tsv"


# ── Aggregate target list for `rule all` ─────────────────────────────────────
def get_all_targets(wildcards):
    """Full default target list (per-sample preprocessing + optional niches +
    annotation-dependent downstream + subcompartments). References producing rules
    by name; resolved at DAG-build time, when every rule is defined."""
    targets = []
    targets += expand(rules.preprocess_umap.output.adata, sample=SAMPLE_IDS)
    targets += expand(rules.preprocess_umap.output.metadata, sample=SAMPLE_IDS)
    targets += expand(rules.preprocess_umap.output.report, sample=SAMPLE_IDS)
    if MODE == "visiumhd":   # HEAD output; only produced when that head runs
        targets += expand(rules.generate_qupath_image.output.qupath_image, sample=SAMPLE_IDS)
    targets += expand(rules.leiden_analysis.output.res_dir,
                      sample=SAMPLE_IDS, resolution=RESOLUTIONS)
    if SPATIAL_NICHES_ENABLED:
        targets.append(rules.spatial_niches.output.concatenated)
        targets.append(rules.spatial_niches.output.plots_dir)
        targets += list(rules.spatial_niches.output.niche_tsvs)
    if HAS_CLUSTER_ANNOT:
        if HAS_INGEST_REF:
            targets += expand(rules.ingest_ref.output.adata_ingested, sample=SAMPLE_IDS)
        targets += expand(rules.annotate_cells.output.adata_annot, sample=SAMPLE_IDS)
        targets += expand(rules.neighbourhood_analysis.output.results_dir,
                          sample=SAMPLE_IDS, annot_type=ANNOT_TYPES)
        targets += expand(rules.pseudobulk_aggregate.output.agg_dir,
                          annot_type=ANNOT_TYPES, analysis_level=ANALYSIS_LEVELS)
        if RUN_DE:
            targets += expand(rules.pseudobulk_de.output.results_dir,
                              annot_type=ANNOT_TYPES, analysis_level=ANALYSIS_LEVELS)
        targets.append(rules.integrate_samples.output.concatenated)
        targets.append(rules.integrate_samples.output.harmony)
        targets.append(rules.integrate_samples.output.sketched)
        targets.append(rules.sample_report.output.report)
        if HAS_GENE_EXPLORATION:
            targets.append(rules.explore_genes_integrated.output.ranges)
            targets += expand(rules.explore_genes_sample.output.done, sample=SAMPLE_IDS)
    if SUBCOMPARTMENTS:
        targets += expand(rules.subcluster.output.sub_dir, subcompartment=SUBCOMPARTMENTS)
    return targets
