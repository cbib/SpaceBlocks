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


def mem_mb_attempt(rule_name, cap_factor=None):
    """mem_mb that grows with the Snakemake retry attempt (base * attempt), so a job
    killed for OOM is automatically resubmitted with more memory. Pairs with the
    global `retries` set in the execution profile. attempt starts at 1 (= base)."""
    base = get_resource(rule_name, "mem_mb")
    def _mem(wildcards, attempt):
        factor = attempt if cap_factor is None else min(attempt, cap_factor)
        return int(base * factor)
    return _mem


# ── Small derivations called from the Snakefile globals block ────────────────
def check_external_annotation():
    """Fail fast (at parse time, before any job) when external_annotation is enabled
    but not fully in place: the annotation column must be present in the metadata for
    EVERY sample. External labels can only come from precomputed_metadata_dir (the
    pipeline's own metadata does not carry them), so that dir is required."""
    cfg = config.get("external_annotation", {}) or {}
    if not cfg.get("enabled", False):
        return
    col = cfg.get("column", "")
    if not col:
        sys.exit("[config error] external_annotation.enabled is true but 'column' is empty.")
    meta_dir = config.get("precomputed_metadata_dir", "") or ""
    if not meta_dir:
        sys.exit("[config error] external_annotation.enabled requires 'precomputed_metadata_dir' "
                 f"to point at a directory of metadata_{{sample}}.tsv files carrying the '{col}' "
                 "column (the pipeline's own metadata does not contain external labels).")
    missing_file, missing_col = [], []
    for s in SAMPLE_IDS:
        f = os.path.join(meta_dir, f"metadata_{s}.tsv")
        if not os.path.isfile(f):
            missing_file.append(f)
            continue
        try:
            header = pd.read_csv(f, sep="\t", nrows=0, comment="#").columns
        except Exception as e:
            sys.exit(f"[config error] could not read external metadata {f}: {e}")
        if col not in header:
            missing_col.append(f)
    if missing_file or missing_col:
        msg = [f"[config error] external_annotation enabled (column '{col}') but not in place:"]
        if missing_file:
            msg.append("  missing metadata file(s):\n    " + "\n    ".join(missing_file))
        if missing_col:
            msg.append(f"  column '{col}' absent in:\n    " + "\n    ".join(missing_col))
        sys.exit("\n".join(msg))


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


# ── HEAD (Xenium 5K) input helper ────────────────────────────────────────────
def xenium_dir_for(sample):
    """Xenium output bundle directory for a sample, from the xenium5k.xenium_dir
    {sample} pattern (e.g. 'data/xenium/{sample}')."""
    pat = (config.get("xenium5k", {}) or {}).get("xenium_dir", "")
    if not pat:
        sys.exit("[config error] mode 'xenium5k' requires xenium5k.xenium_dir "
                 "(a {sample} pattern to each Xenium output bundle).")
    return pat.format(sample=sample)

# ── HEAD (Atera) ─────────────────────────────────────────────────────────────
def atera_dir_for(sample):
    """Atera output bundle directory for a sample, from the atera.atera_dir
    {sample} pattern (e.g. 'data/atera/{sample}/outs')."""
    pat = (config.get("atera", {}) or {}).get("atera_dir", "")
    if not pat:
        sys.exit("[config error] mode 'atera' requires atera.atera_dir "
                 "(a {sample} pattern to each Atera outs/ bundle).")
    return pat.format(sample=sample)


def he_file_for(sample, key, required=True):
    """Resolve one of the OPTIONAL H&E supplemental files (he_image, he_alignment,
    he_keypoints) for a sample. These ship as separate downloads from the outs/
    bundle, so the config gives absolute {sample} patterns rather than the head
    discovering them. Returns "" when unset and required=False — used for
    he_keypoints, whose absence only disables the alignment QA."""
    pat = (config.get("atera", {}) or {}).get(key, "")
    if not pat:
        if required:
            sys.exit(f"[config error] atera.{key} is required when the optional H&E "
                     "annotation image is enabled (atera.he_image + atera.he_alignment).")
        return ""
    return pat.format(sample=sample)


def _ate_prepare_inputs(wildcards):
    """Inputs for prepare_input_ate. The H&E scale-factor JSON and background image are
    required only when the optional H&E QuPath rule is active, so they are added
    conditionally — a static input: block cannot express that."""
    inputs = {
        "done": rules.convert_zarr_ate.output.done.format(sample=wildcards.sample),
        "qupath_meta": rules.generate_qupath_ate.output.qupath_meta.format(
            sample=wildcards.sample),
    }
    if ATERA_HAS_HE:
        # Only the embedded background is consumed here; the raw-H&E QuPath image and its
        # polygon affine were dropped, so he_meta is no longer an input.
        inputs["he_background"] = rules.generate_qupath_he_ate.output.he_background.format(
            sample=wildcards.sample)
    # Track the region GeoJSON when present, so editing/renaming it retriggers the
    # contract build. It is read by filename (regions are optional), not hard-required —
    # same "declare only when present" pattern as _preprocess_inputs' precomputed_meta.
    for _suffix in ("_he_background.geojson", "_morphology.geojson", "_he.geojson"):
        _gj = os.path.join(GEOJ_DIR, f"{wildcards.sample}{_suffix}")
        if os.path.isfile(_gj):
            inputs["geojson"] = _gj
            break
    return inputs

# ── HEAD (MERSCOPE) input helper ─────────────────────────────────────────────
def merscope_dir_for(sample):
    """MERSCOPE region directory for a sample, from the merscope.merscope_dir
    {sample} pattern (e.g. 'data/merscope/{sample}'). The region dir holds the two
    Vizgen CSVs (cell_by_gene, cell_metadata) and the images/ folder (mosaic TIFFs +
    micron_to_mosaic_pixel_transform.csv)."""
    pat = (config.get("merscope", {}) or {}).get("merscope_dir", "")
    if not pat:
        sys.exit("[config error] mode 'merscope' requires merscope.merscope_dir "
                 "(a {sample} pattern to each MERSCOPE region directory).")
    return pat.format(sample=sample)


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
        "validation": rules.validate_input.output.report.format(sample=wc.sample),
    }
    if USE_PRECOMPUTED:
        precomp_dir = config.get("precomputed_metadata_dir", "")
        if precomp_dir:
            meta_file = os.path.join(precomp_dir, f"metadata_{wc.sample}.tsv")
            if os.path.isfile(meta_file):
                inputs["precomputed_meta"] = meta_file
    th_tsv = config.get("per_sample_qc", "") or ""
    if th_tsv and os.path.isfile(th_tsv):
        inputs["thresholds_tsv"] = th_tsv
    # External-annotation mask (keep_unannotated=false): the external metadata drives
    # the cell set, so it is a tracked input of preprocessing.
    _ext = config.get("external_annotation", {}) or {}
    if _ext.get("enabled") and not _ext.get("keep_unannotated", True):
        _meta_dir = config.get("precomputed_metadata_dir", "") or ""
        if _meta_dir:
            inputs["external_meta"] = os.path.join(_meta_dir, f"metadata_{wc.sample}.tsv")
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
        return rules.ingest_ref.output.adata_ingested.format(sample=wc.sample)
    return rules.preprocess_umap.output.adata.format(sample=wc.sample)


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
    targets += list(QUPATH_IMAGES)   # head QuPath image(s) for the active mode (or [])
    if RUN_LEIDEN_ANALYSIS:
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
