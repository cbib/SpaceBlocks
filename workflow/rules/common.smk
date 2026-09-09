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


# ── Expand {base_dir} in config paths ────────────────────────────────────────
BASE_DIR = config.get("base_dir", "")


def _expand_base_dir(node, root="config"):
    """
    Expand config paths (recursively), replacing only {base_dir}
    leaving other placeholders (for example, the wildcard {sample}) intact.
    """
    if isinstance(node, str):
        if "{base_dir}" in node:
            expanded = node.replace("{base_dir}", BASE_DIR)
            print(f"[config] expanded {root}: {expanded}", file=sys.stderr)
            return expanded
        return node
    if isinstance(node, dict):
        return {k: _expand_base_dir(v, f"{root}.{k}") for k, v in node.items()}
    if isinstance(node, list):
        return [_expand_base_dir(v, f"{root}[{i}]") for i, v in enumerate(node)]
    return node


# ── Small derivations called from the Snakefile globals block ────────────────
def check_external_annotation():
    """Fail fast (at parse time, before any job) when external_annotation is enabled
    but not fully in place. Metadata TSV columns are checked here; in decoupled mode,
    an absent metadata directory means the configured column will instead be checked
    in each contract h5ad by validate_input."""
    cfg = config.get("external_annotation", {}) or {}
    if not cfg.get("enabled", False):
        return
    col = cfg.get("column", "")
    if not col:
        sys.exit(
            "[config error] external_annotation.enabled is true but 'column' is empty."
        )
    meta_dir = config.get("precomputed_metadata_dir", "") or ""
    if not meta_dir:
        if IS_DECOUPLED:
            return
        sys.exit(
            "[config error] external_annotation.enabled requires 'precomputed_metadata_dir' "
            f"to point at a directory of metadata_{{sample}}.tsv files carrying the '{col}' "
            "column. In decoupled mode only, the column may instead be stored directly "
            "in each contract h5ad's obs."
        )
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
        msg = [
            f"[config error] external_annotation enabled (column '{col}') but not in place:"
        ]
        if missing_file:
            msg.append(
                "  missing metadata file(s):\n    " + "\n    ".join(missing_file)
            )
        if missing_col:
            msg.append(
                f"  column '{col}' absent in:\n    " + "\n    ".join(missing_col)
            )
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
            cols.append(c)
            vals.append(str(v))
            colors.append(str(hexc))
    return cols, vals, colors


def _validate_pseudobulk_analysis(spec):
    """Validate pseudobulk settings that JSON Schema cannot check across fields."""
    if not isinstance(spec, dict):
        sys.exit(
            "[config error] each analysis.pseudobulk.analyses entry must be a mapping."
        )

    name = str(spec.get("name", ""))
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        sys.exit(
            "[config error] pseudobulk analysis names may contain only letters, "
            f"numbers, '_' and '-'; got '{name}'."
        )

    aggregation = str(spec.get("aggregation", ""))
    if aggregation not in _PSEUDOBULK_AGGREGATIONS:
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' has unsupported aggregation "
            f"'{aggregation}'. Choose one of {sorted(_PSEUDOBULK_AGGREGATIONS)}."
        )

    group_by = list(spec.get("group_by") or [])
    if not group_by or any(not isinstance(c, str) or not c for c in group_by):
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' needs at least one "
            "non-empty group_by column."
        )
    if len(group_by) != len(set(group_by)):
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' repeats a group_by column."
        )

    paired_by = spec.get("paired_by", "") or ""
    if not isinstance(paired_by, str):
        sys.exit(
            f"[config error] pseudobulk analysis '{name}'.paired_by must be a string."
        )

    covariates = list(spec.get("covariates") or [])
    if any(not isinstance(c, str) or not c for c in covariates):
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' has an invalid covariate."
        )
    if len(covariates) != len(set(covariates)):
        sys.exit(f"[config error] pseudobulk analysis '{name}' repeats a covariate.")
    overlap = sorted(set(group_by) & set(covariates))
    if overlap:
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' lists {overlap} as both "
            "group_by and covariates."
        )
    if paired_by and paired_by in group_by:
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' lists '{paired_by}' as both "
            "paired_by and group_by."
        )
    if paired_by and paired_by in covariates:
        sys.exit(
            f"[config error] pseudobulk analysis '{name}' lists '{paired_by}' as both "
            "paired_by and a covariate."
        )

    exclude_levels = spec.get("exclude_levels") or {}
    if not isinstance(exclude_levels, dict):
        sys.exit(
            f"[config error] pseudobulk analysis '{name}'.exclude_levels must be a mapping."
        )

    contrast_names = set()
    for contrast in list(spec.get("contrasts") or []):
        if not isinstance(contrast, dict):
            sys.exit(
                f"[config error] contrasts in pseudobulk analysis '{name}' must be mappings."
            )
        contrast_name = str(contrast.get("name", ""))
        if not re.fullmatch(r"[A-Za-z0-9_-]+", contrast_name):
            sys.exit(
                f"[config error] contrast names in pseudobulk analysis '{name}' may "
                "contain only letters, numbers, '_' and '-'."
            )
        if contrast_name in contrast_names:
            sys.exit(
                f"[config error] duplicate contrast '{contrast_name}' in '{name}'."
            )
        contrast_names.add(contrast_name)
        for side in ("numerator", "denominator"):
            values = contrast.get(side)
            if not isinstance(values, dict) or set(values) != set(group_by):
                sys.exit(
                    f"[config error] contrast '{contrast_name}' {side} must specify "
                    f"exactly the group_by columns {group_by}."
                )
            excluded_values = {
                column: {str(value) for value in exclude_levels.get(column, [])}
                for column in group_by
            }
            conflicts = [
                f"{column}={value}"
                for column, value in values.items()
                if str(value) in excluded_values[column]
            ]
            if conflicts:
                sys.exit(
                    f"[config error] contrast '{contrast_name}' {side} in '{name}' "
                    f"uses excluded level(s): {conflicts}."
                )
        if contrast["numerator"] == contrast["denominator"]:
            sys.exit(
                f"[config error] contrast '{contrast_name}' in '{name}' has identical "
                "numerator and denominator groups."
            )

    return {
        **spec,
        "name": name,
        "aggregation": aggregation,
        "group_by": group_by,
        "paired_by": paired_by,
        "covariates": covariates,
        "exclude_levels": dict(exclude_levels),
        "contrasts": list(spec.get("contrasts") or []),
        "lrt": dict(spec.get("lrt") or {}),
    }


def _legacy_pseudobulk_analyses():
    """Translate the former region-only config for existing demo/user configs.

    New configurations should use analysis.pseudobulk.analyses. The unimplemented
    by_niche_region value is deliberately rejected rather than carried forward.
    """
    levels = list(ANALYSIS.get("analysis_levels", ["by_region", "by_celltype_region"]))
    unsupported = sorted(set(levels) - {"by_region", "by_celltype_region"})
    if unsupported:
        sys.exit(
            "[config error] unsupported legacy pseudobulk analysis_levels: "
            f"{unsupported}. 'by_niche_region' has been removed because it was never implemented."
        )
    contrasts = [
        {
            "name": re.sub(r"[^A-Za-z0-9_-]+", "_", f"{a}_vs_{b}").strip("_"),
            "numerator": {"region_annotation": str(a)},
            "denominator": {"region_annotation": str(b)},
        }
        for a, b in combinations(REGION_LEVELS, 2)
    ]
    return [
        {
            "name": level,
            "aggregation": level,
            "group_by": ["region_annotation"],
            "paired_by": "sample",
            "covariates": [],
            "exclude_levels": {"region_annotation": ["Unlabeled", "Bubble"]},
            "contrasts": contrasts,
            "lrt": {"enabled": len(REGION_LEVELS) > 2},
        }
        for level in levels
    ]


def pseudobulk_analysis(name):
    """Return one validated named pseudobulk analysis from the Snakefile globals."""
    return PSEUDOBULK_ANALYSIS_BY_NAME[str(name)]


def pseudobulk_group_label(spec, values):
    """Stable display/factor label for one configured group combination."""
    columns = spec["group_by"]
    if len(columns) == 1:
        return str(values[columns[0]])
    return " | ".join(f"{column}={values[column]}" for column in columns)


def pseudobulk_contrast_param(analysis_name, field):
    """Flatten named contrast mappings for safe Snakemake Python-to-R transfer."""
    spec = pseudobulk_analysis(analysis_name)
    if field == "names":
        return [str(c["name"]) for c in spec["contrasts"]]
    if field not in {"numerator", "denominator"}:
        raise ValueError(f"Unknown pseudobulk contrast field: {field}")
    return [pseudobulk_group_label(spec, c[field]) for c in spec["contrasts"]]


def pseudobulk_group_palette(analysis_name):
    """Palette for a one-column comparison; combined groups use R fallbacks."""
    spec = pseudobulk_analysis(analysis_name)
    if len(spec["group_by"]) != 1:
        return {}
    column = spec["group_by"][0]
    if column == "region_annotation":
        return ANALYSIS.get("region_colors", {}) or {}
    return SAMPLE_COLORS.get(column, {}) if isinstance(SAMPLE_COLORS, dict) else {}


def pseudobulk_condition_order(analysis_name):
    """Preferred heatmap order for a region-only comparison."""
    spec = pseudobulk_analysis(analysis_name)
    if spec["group_by"] == ["region_annotation"]:
        return [str(level) for level in REGION_LEVELS]
    return []


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
        sys.exit(
            "[config error] mode 'xenium5k' requires xenium5k.xenium_dir "
            "(a {sample} pattern to each Xenium output bundle)."
        )
    return pat.format(sample=sample)


# ── HEAD (Atera) ─────────────────────────────────────────────────────────────
def atera_dir_for(sample):
    """Atera output bundle directory for a sample, from the atera.atera_dir
    {sample} pattern (e.g. 'data/atera/{sample}/outs')."""
    pat = (config.get("atera", {}) or {}).get("atera_dir", "")
    if not pat:
        sys.exit(
            "[config error] mode 'atera' requires atera.atera_dir "
            "(a {sample} pattern to each Atera outs/ bundle)."
        )
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
            sys.exit(
                f"[config error] atera.{key} is required when the optional H&E "
                "annotation image is enabled (atera.he_image + atera.he_alignment)."
            )
        return ""
    return pat.format(sample=sample)


def _ate_prepare_inputs(wildcards):
    """Inputs for prepare_input_ate. The H&E scale-factor JSON and background image are
    required only when the optional H&E QuPath rule is active, so they are added
    conditionally — a static input: block cannot express that."""
    inputs = {
        "done": rules.convert_zarr_ate.output.done.format(sample=wildcards.sample),
        "qupath_meta": rules.generate_qupath_ate.output.qupath_meta.format(
            sample=wildcards.sample
        ),
    }
    if ATERA_HAS_HE:
        # Only the embedded background is consumed here; the raw-H&E QuPath image and its
        # polygon affine were dropped, so he_meta is no longer an input.
        inputs["he_background"] = (
            rules.generate_qupath_he_ate.output.he_background.format(
                sample=wildcards.sample
            )
        )
    # Region annotations are optional. When one is present, track both the source
    # and its validation report so the contract is rebuilt only after validation.
    _gj = _find_geojson(wildcards.sample)
    if _gj:
        inputs["geojson"] = _gj
        inputs["geojson_validation"] = _geojson_validation_report(wildcards.sample)
    return inputs


def _geojson_suffixes():
    """Accepted region-annotation filenames for the active head, in priority order."""
    return {
        "visiumhd": ("_tissue_hires_image.geojson",),
        "xenium5k": ("_morphology.geojson",),
        "atera": ("_he_background.geojson", "_morphology.geojson"),
        "merscope": ("_morphology.geojson",),
    }.get(MODE, ())


def _find_geojson(sample):
    """Resolve an optional region GeoJSON using filenames supported by the active head."""
    if not GEOJ_DIR or not _geojson_suffixes():
        return None

    accepted = []
    for _suffix in _geojson_suffixes():
        accepted.append(os.path.join(GEOJ_DIR, f"{sample}{_suffix}"))

    # A sample-prefixed GeoJSON with any other suffix is almost certainly a typo or
    # an export for a different coordinate system. Fail before silently ignoring it.
    discovered = sorted(glob.glob(os.path.join(GEOJ_DIR, f"{sample}_*.geojson")))
    unsupported = [path for path in discovered if path not in accepted]
    if unsupported:
        sys.exit(
            f"[config error] Unsupported GeoJSON filename(s) for sample '{sample}' "
            f"in mode '{MODE}': {unsupported}. Expected one of: {accepted}"
        )
    return next((path for path in accepted if os.path.isfile(path)), None)


def _geojson_validation_report(sample):
    """Validation marker consumed by a head only when that sample has a GeoJSON."""
    if _find_geojson(sample) is None:
        return []
    return f"{SAMPLES_DIR}/{sample}/validation/geojson_validation.json"


def _external_metadata_input(wildcards):
    """Tracked per-sample external metadata, or no input when the feature is off."""
    cfg = config.get("external_annotation", {}) or {}
    metadata_dir = config.get("precomputed_metadata_dir", "") or ""
    if not cfg.get("enabled", False) or not metadata_dir:
        return []
    return os.path.join(
        metadata_dir,
        f"metadata_{wildcards.sample}.tsv",
    )


# ── HEAD (MERSCOPE) input helper ─────────────────────────────────────────────
def merscope_dir_for(sample):
    """MERSCOPE region directory for a sample, from the merscope.merscope_dir
    {sample} pattern (e.g. 'data/merscope/{sample}'). The region dir holds the two
    Vizgen CSVs (cell_by_gene, cell_metadata) and the images/ folder (mosaic TIFFs +
    micron_to_mosaic_pixel_transform.csv)."""
    pat = (config.get("merscope", {}) or {}).get("merscope_dir", "")
    if not pat:
        sys.exit(
            "[config error] mode 'merscope' requires merscope.merscope_dir "
            "(a {sample} pattern to each MERSCOPE region directory)."
        )
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
            inputs["external_meta"] = os.path.join(
                _meta_dir, f"metadata_{wc.sample}.tsv"
            )
    return inputs


def _spatial_niches_inputs(wc):
    """Inputs for spatial_niches. When use_precomputed is set with an external
    niche_dir, the per-sample niche_{sample}.tsv are added as MANDATORY inputs so
    Snakemake aborts before the job if any is absent (complementary to the in-script
    safeguard, which also catches empty/partial TSVs). When niche_dir is empty the
    rule reloads from its own output dir, so nothing is added here (no circular
    dependency)."""
    inputs = {
        "adatas": expand(
            f"{SAMPLES_DIR}/{{sample}}/adata_{{sample}}.h5ad", sample=SAMPLE_IDS
        ),
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
    targets += expand(rules.process_geojson.output.report, sample=SAMPLES_WITH_GEOJSON)
    targets += list(QUPATH_IMAGES)  # head QuPath image(s) for the active mode (or [])
    if RUN_LEIDEN_ANALYSIS:
        targets += expand(
            rules.leiden_analysis.output.res_dir,
            sample=SAMPLE_IDS,
            resolution=RESOLUTIONS,
        )
    if SPATIAL_NICHES_ENABLED:
        targets.append(rules.spatial_niches.output.concatenated)
        targets.append(rules.spatial_niches.output.plots_dir)
        targets += list(rules.spatial_niches.output.niche_tsvs)
    if HAS_CLUSTER_ANNOT:
        if HAS_INGEST_REF:
            targets += expand(rules.ingest_ref.output.adata_ingested, sample=SAMPLE_IDS)
        targets += expand(rules.annotate_cells.output.adata_annot, sample=SAMPLE_IDS)
        targets += expand(
            rules.neighbourhood_analysis.output.results_dir,
            sample=SAMPLE_IDS,
            annot_type=ANNOT_TYPES,
        )
        targets += expand(
            rules.pseudobulk_aggregate.output.agg_dir,
            annot_type=ANNOT_TYPES,
            analysis_name=PSEUDOBULK_ANALYSIS_NAMES,
        )
        if RUN_DE:
            targets += expand(
                rules.pseudobulk_de.output.results_dir,
                annot_type=ANNOT_TYPES,
                analysis_name=PSEUDOBULK_ANALYSIS_NAMES,
            )
        targets.append(rules.integrate_samples.output.concatenated)
        targets.append(rules.integrate_samples.output.harmony)
        targets.append(rules.integrate_samples.output.sketched)
        targets.append(rules.sample_report.output.report)
        if HAS_GENE_EXPLORATION:
            targets.append(rules.explore_genes_integrated.output.ranges)
            targets += expand(rules.explore_genes_sample.output.done, sample=SAMPLE_IDS)
    if SUBCOMPARTMENTS:
        targets += expand(
            rules.subcluster.output.sub_dir, subcompartment=SUBCOMPARTMENTS
        )
    return targets
