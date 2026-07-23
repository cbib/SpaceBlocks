# SpaceBlocks — pipeline documentation

A modular Snakemake pipeline for **single-cell-resolution spatial transcriptomics**.
A technology-specific **Headblock** turns raw platform output into one standardized
*contract* file; a technology-agnostic **Coreblock** runs the entire analysis on that
contract, in three steps (preprocessing → postprocessing → exploration).

> **Scope.** The goal is a Coreblock reusable across single-cell-resolution spatial
> platforms (Visium HD today; others by adding a Headblock). It is **not** designed to
> integrate or compare data *across* technologies — each run analyses one platform.

---

## 1. Architecture

```
── Headblock (Visium HD specific) ───────────────────────────────────────────────
spaceranger_count                 Space Ranger 4.0.1 (installed by the user)
      │
generate_qupath_image             copy the hires tissue image out of the SR tree
      │   ‖ MANUAL BREAK: annotate regions in QuPath, export
      ‖     {sample}_tissue_hires_image.geojson → geojson_path
      │
prepare_input                     raw 10x matrix + cell centroids + embedded image
      ▼   writes the UNFILTERED CONTRACT h5ad  (X = raw counts, no filtering)
── Coreblock (technology-agnostic) ──────────────────────────────────────────────
STEP 1 · preprocessing
  validate_input → [qc_sweep] → preprocess_umap → leiden_analysis
                                → generate_annotation_template
                                → ingest_ref (optional) → spatial_niches (optional)
STEP 2 · postprocessing
  annotate_cells → integrate_samples → pseudobulk_aggregate → pseudobulk_de
                 → neighbourhood_analysis → subcluster → sample_report
STEP 3 · exploration
  explore_genes_integrated → explore_genes_sample
```

The Headblock is the only technology-specific part. **Adding a new platform = write a new
Headblock** that emits the contract h5ad; the Coreblock is untouched. See `docs/rules.md` for a
one-paragraph description of every rule.

---

## 2. The contract (Headblock → Coreblock hand-off)

The single coupling point between Headblock and Coreblock is one h5ad per sample.

- **Path** = `config["contract"]["unfiltered_h5ad"]`, derived in the Snakefile from
  `SAMPLES_DIR` as `{SAMPLES_DIR}/{sample}/{sample}_unfiltered.h5ad`. `prepare_input`
  WRITES it; `validate_input`, `qc_sweep`, and `preprocess_umap` READ it via the same
  config key, so producer and consumers can never drift.
- **Shape** (written by `prepare_input`):
  - `X` — raw integer counts (no filtering, no normalisation).
  - `obsm["spatial"]` — per-cell `(x, y)` centroids, finite, shape `(n, ≥2)`.
  - `uns["spatial"][sample]` — `{images: {hires}, scalefactors}` when a tissue image
    exists; omitted for image-free platforms (the Coreblock then scatters on coordinates).
  - `obs` — `sample`, `cell_id`, and `region_annotation` when a GeoJSON was provided.
- **The image is a Headblock concern.** The Coreblock never rebuilds an image; it consumes the
  embedded `uns["spatial"]` image when present and otherwise scatters on `obsm["spatial"]`
  with an explicit `spot_size` (cell centroids are far smaller than the scalefactor-derived
  Visium spot size, so `spot_size` is always passed).

`config["outdir"]` is set to `SAMPLES_DIR` and is the per-sample output root for
`validate_input` (`{outdir}/{sample}/validation/`) and `qc_sweep` (`{outdir}/{sample}/qc_sweep/`).

---

## 3. Requirements & installation

- **Snakemake ≥ 8** and **conda/mamba** — one conda env per rule group (`--use-conda`).
- **SLURM executor plugin** — Snakemake 8 moved cluster submission to a plugin:
  `pip install snakemake-executor-plugin-slurm` (the profile sets `executor: slurm`).
- **Space Ranger 4.0.1** — installed and run by the user (proprietary; path in
  `config["spaceranger"]`). It cannot be conda/container-packaged.
- **QuPath** — for the manual region-annotation step (exports GeoJSON).

```bash
git clone <GITHUB_URL> && cd <repo>
mamba create -n spaceblocks -c conda-forge -c bioconda snakemake mamba
mamba activate spaceblocks
pip install snakemake-executor-plugin-slurm
```

---

## 4. Inputs

Under `config/`:

- **`samples.csv`** — Headblock sheet (`sample_id`, `fastq_dir`, `slide`, `area`, and optional
  re-sequencing dirs). Validated against `schemas/samples.schema.yaml`.
- **`core_samples.tsv`** — the Coreblock anchor sheet. Required column `sample`; optional
  metadata columns (e.g. `patient`, `type`) that are stamped into `obs` and surfaced
  downstream. This sheet drives the Coreblock sample list; a `sample` value need not equal
  `obs["sample"]` inside the contract (samples may be renamed for publication).
- **Marker / annotation / query TSVs** — `snakemake_cell_markers.tsv`,
  `annotation_markers.tsv`, `cluster_annotations.tsv`, `gene_queries.tsv`.
- **QuPath GeoJSON** — one `{sample}_tissue_hires_image.geojson` per sample in
  `geojson_path`, exported after the manual annotation step.
- **`config.yaml`** — all parameters (validated against `schemas/config.schema.yaml`).

---

## 5. Configuration highlights

`config/config.yaml` is the source of truth; the key blocks are:

```yaml
mode: "visiumhd"                    # a Headblock technology name, or "decoupled"

contract:                          # semantic keys only (path/outdir set in Snakefile)
  sample_key: "sample"
  spatial_key: "spatial"
  require_region: false
  require_raw_counts: true
  mito_prefix: ["MT-", "mt-"]

core_samples: "config/core_samples.tsv"
contract_dir: ""                   # decoupled mode: where to READ pre-existing contracts

ingest_ref: "resources/reference.h5ad"   # optional scanpy label transfer (shared by
ingest_ref_label_key: "cell_type"        #   ingest_ref and the qc_sweep overlay)

integration:
  integrate_key: "sample"          # obs column Harmony corrects on
extra_annotations:
  columns: ["patient", "type"]     # core-sheet columns surfaced downstream (NOT DE covariates)
sample_colors:                     # column -> {value -> hex}; grey (#cccccc) fallback
  patient: { "Patient 1": "#000000", ... }

analysis:
  resolution_scan_min: 0.8         # Leiden resolution ladder
  resolution_scan_max: 1.2
  resolution_scan_step: 0.2
  run_pseudobulk_de: false         # DESeq2 step is optional
  region_levels: [...]             # ordered region names (must match the GeoJSON values)
  region_colors: {...}

qc_sweep:                          # OPTIONAL diagnostic (target: qc_sweep_all)
  ingest_enabled: false            # overlay cell types using the top-level ingest_ref
  dpi: 300
  thresholds:                      # LISTS of candidate cut-offs to sweep
    min_genes: [50, 100, 200]
    min_counts: [1, 10, 50]
    max_counts: [3000, 5000]
    max_pct_mt: [10, 15, 20]

preprocess_thresholds: ""          # OPTIONAL per-sample QC TSV (overrides analysis.* per sample)

resources:                         # per-rule mem_mb / runtime / threads
  prepare_input: { mem_mb: 32768, runtime: 240, threads: 4 }
  qc_sweep:      { mem_mb: 49152, runtime: 180, threads: 4 }
  # ...
```

**Per-sample QC thresholds.** `preprocess_thresholds` may point to a TSV (first column =
sample, columns = `min_counts`/`min_genes`/`min_cells`/optional `max_counts`/`max_pct_mt`);
values override the `analysis.*` defaults key-by-key, with absent samples/columns falling
back to config. This override lives in `preprocess_umap`, where filtering happens — NOT in
`qc_sweep`, whose thresholds are global candidate lists.

**Schema note.** `config.schema.yaml` validates config BEFORE the Snakefile injects
`contract.unfiltered_h5ad`/`outdir`, so those two are intentionally absent from the schema.
The `contract`, `qc_sweep`, `integration`, `extra_annotations`, `sample_colors`, and
`resources` keys ARE validated.

**External annotation.** Cell-type labels produced by an external tool can be overlaid via
`external_annotation` (`enabled` + `column`) plus a `precomputed_metadata_dir` of
`metadata_{sample}.tsv` files carrying that column (barcode-indexed). `validate_input`
checks barcode overlap up front, and `annotate_cells` writes `obs["cell_type_external"]`.
Two behaviours matter:

- **It takes over.** When enabled, the external annotation becomes the *primary* annotation
  everywhere — `annotation_types` (neighbourhood/pseudobulk) and the colour column for
  `explore_genes`/`subcluster` both default to it, and `tsv_annotation` is **not** analysed
  by default. This is deliberate: otherwise cells the external tool left unannotated would be
  re-labelled by the tsv cluster mapping and leak as noise into downstream analyses (e.g.
  subcompartments). To compare hand-made vs external, set `annotation_types:
  ["tsv_annotation", "external_annotation"]` explicitly.
- **`keep_unannotated` decides QC.** `true` (default) keeps normal QC and overlays labels
  (unmatched cells → `Unannotated`); `false` keeps *only* the externally-annotated cells and
  **skips the pipeline QC entirely** — the external labels are the QC decision — so niches,
  annotation, and all downstream run on exactly that cell set.

---

## 6. Running

The workflow has a **manual break** for region annotation. Run it in stages; always
dry-run (`-n`) first and pass the profile (`--profile profiles/default`).

```bash
# 0. Configure config.yaml, samples.csv, core_samples.tsv.

# 1. Headblock: build the hires images to annotate.
snakemake qupath_images --profile profiles/default
#    → open each Samples/<sample>/QuPath_image/*.png in QuPath, draw regions,
#      export {sample}_tissue_hires_image.geojson into geojson_path.

# 2. (Optional) sweep QC thresholds before committing filters.
snakemake qc_sweep_all --profile profiles/default

# 3. STEP 1 — preprocessing.
snakemake run_preprocessing --profile profiles/default

# 4. STEP 2 — postprocessing (after filling cluster_annotations).
snakemake run_postprocessing --profile profiles/default

# 5. STEP 3 — exploration.
snakemake run_exploration --profile profiles/default
```

**Convenience targets:** `qupath_images`, `qc_sweep_all`, `run_preprocessing`,
`run_postprocessing`, `run_exploration`, `subcluster_all`, `explore_genes`, and `all`
(everything). Target input lists reference producing rules symbolically
(`rules.<rule>.output.<name>`), so a path change in a rule propagates automatically.

**Run modes.**
- `mode: visiumhd` — the Visium HD Headblock builds the contracts in this run.
- `mode: decoupled` — the Coreblock consumes PRE-EXISTING contracts from `contract_dir`
  (built by any Headblock, possibly in a separate run); no Headblock rules are included and
  validation aborts up front if any sample is missing its contract.

**Execution & robustness (SLURM profile).** `profiles/default/config.yaml` sets
`executor: slurm`, `jobs: 30`, `use-conda: true`, `keep-going: true`,
`rerun-incomplete: true`, `latency-wait: 60`, and `retries: 3`. On a failure a job is
resubmitted up to three times, and **memory grows with the attempt** (`mem_mb = base ×
attempt`) so an OOM-killed job is retried with more RAM instead of failing identically.
Because memory escalates, make sure `base_mem × max_attempts` for your heaviest rules
still fits your nodes (lower `retries` or the base `resources.*.mem_mb` otherwise);
runtime is not scaled, to avoid exceeding walltime partitions. No rule runs on the head
node (no `localrule`s), so head-node-restricted clusters are supported.

---

## 7. Output structure

```
{post_processing_outdir}/
├── Samples/{sample}/
│   ├── {sample}_unfiltered.h5ad                 CONTRACT (prepare_input)
│   ├── QuPath_image/…_tissue_hires_image.png    generate_qupath_image
│   ├── validation/input_validation.json         validate_input
│   ├── qc_sweep/                                 qc_sweep (optional)
│   ├── adata_{sample}.h5ad, metadata_{sample}.tsv, {sample}_report.tsv   preprocess_umap
│   ├── leiden_resolution_*/                      leiden_analysis
│   ├── annotation/, adata_{sample}_annotated.h5ad   annotate_cells
│   └── neighbourhood_analysis/{annot_type}/      neighbourhood_analysis
├── integrated_samples/                           concatenated / harmony / sketched / samples_report.pdf
├── pseudobulk/{annot_type}/{analysis_level}/     aggregated/ + de_results/
├── Subcompartments/{name}/{Harmony,NoHarmony}/   subcluster
├── spatial_niches/                               tsv/, plots/ (+ individual_niches/)
├── gene_exploration/                             expression_ranges.tsv, {entry}/…
└── cluster_annotations_template.tsv              generate_annotation_template
```

---

## 8. Key design decisions

- **Headblock/Coreblock split.** Technology-specific code produces a standardized unfiltered
  contract; the common Coreblock consumes it. New platform = new Headblock, Coreblock unchanged.
- **Contract path is one config key**, referenced by producer and all consumers.
- **`validate_input` is a DAG gate.** A hard failure raises → no JSON report → the core
  cannot start for that sample. Soft issues (missing region annotation, no mito genes,
  no image) are recorded, not fatal.
- **`qc_sweep` is diagnostic only** — never filters, clusters, or writes an h5ad; one
  landing file per parameter; violins rasterise only the per-cell strip; the ingest
  overlay reuses the top-level `ingest_ref`.
- **`generate_qupath_image` separate from `prepare_input`** — a cheap copy meant to be
  grouped with `spaceranger` in a later DAG revision (the GeoJSON must precede
  `prepare_input`).
- **`preprocess_umap` is gated on the validation report** and reads the contract h5ad; it
  is technology-agnostic.
- **Colour palettes are config-driven** (regions, sample metadata) with a grey fallback,
  applied consistently across UMAP/spatial/violin/heatmap outputs.

---

## 9. Extending to a new platform

Write a new Headblock (e.g. `rules/head_xenium5k` rules + a `xenium5k` entry in `KNOWN_HEADS`)
that emits the contract h5ad described in §2, then set `mode: xenium5k`. The entire Coreblock
(steps 1–3) runs unchanged. Alternatively, prepare contracts with any external tool and
run the Coreblock in `mode: decoupled`.
