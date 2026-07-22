# Design & architecture

SpaceBlocks is built around one idea: a technology-specific **HeadBlock** turns raw data into a standardized *contract* object, and a technology-agnostic **CoreBlock** does all the analysis on top of it.

Adding a platform means writing a new HeadBlock; the CoreBlock never changes.

Importantly, the CoreBlocks immediately admit external data, as soon as they pass the `validate_input` rule (see [below](design.md#the-validation-contract)).

## Repository structure

```
SpaceBlocks/
├── workflow/                 the executable workflow (auto-found by Snakemake)
│   ├── Snakefile             globals, mode selection, targets
│   ├── rules/*.smk           one file per rule (+ common.smk for shared helpers)
│   ├── scripts/*.py, *.R     rule implementations
│   ├── envs/*.yaml           one Conda env per rule group
│   └── schemas/*.yaml        config + sample-sheet validation
├── config/                   config.yaml, README.md (config reference), sample sheets
├── docs/                     this documentation site (MkDocs)
├── .tests/                   tiny synthetic decoupled dataset for CI
├── reproduction/             full public-data runs (data fetched, not committed)
├── tools/                    stand-alone helper scripts
└── profiles/default/         SLURM profile (retries, resources)
```

## The Blocks

```
── HeadBlock (technology-specific) ──────────────────────────────────────────────
  mode: visiumhd      spaceranger_count_vhd → generate_qupath_vhd → prepare_input_vhd
  mode: xenium5k      convert_zarr_x5k      → generate_qupath_x5k  → prepare_input_x5k
  mode: decoupled     (no HeadBlock — you provide the contract h5ads)
        │  ‖ MANUAL BREAK: annotate regions in QuPath, export GeoJSON → geojson_path
        ▼  each HeadBlock WRITES the UNFILTERED CONTRACT h5ad (X = raw counts)
── CoreBlock (technology-agnostic) ──────────────────────────────────────────────
  STEP 1 · preprocessing
    validate_input → [qc_sweep] → preprocess_umap → leiden_analysis
                                 → generate_annotation_template
                                 → ingest_ref (opt.) → spatial_niches (opt.)
  STEP 2 · postprocessing
    annotate_cells → integrate_samples → pseudobulk_aggregate → pseudobulk_de
                   → neighbourhood_analysis → subcluster → sample_report
  STEP 3 · exploration
    explore_genes_integrated → explore_genes_sample
```

See the [Rule reference](rules.md) for a one-paragraph description of every rule.

## The validation contract

The point where the HeadBlocks meet the CoreBlocks is the `validate_input` rule. This validation rule writes a `.json` contract that confirms the data is suitable for entering pre- and postprocessing steps.

It is a single coupling point with one h5ad per sample. Every HeadBlock produces the same shape, so the CoreBlock is identical for all technologies:

- **Path** — one config key (`contract.unfiltered_h5ad`, derived in the Snakefile), *written* by
  the active HeadBlock's `prepare_input_*` and *read* by `validate_input`, `qc_sweep`, and
  `preprocess_umap`, so producer and consumers can never drift.
- **Shape:**
  - `X` — raw integer counts (no filtering, no normalisation).
  - `obsm["spatial"]` — per-cell `(x, y)` centroids, finite.
  - `uns["spatial"][sample]` — `{images: {hires}, scalefactors}` when a tissue/morphology image
    exists; omitted otherwise (the CoreBlock then scatters on coordinates).
  - `obs` — `sample`, `cell_id`, and `region_annotation` when a GeoJSON was provided.
- **The image is a HeadBlock concern.** The CoreBlock never rebuilds an image — it consumes the
  embedded `uns["spatial"]` image when present and otherwise scatters on `obsm["spatial"]`.

## Output structure

```
{post_processing_outdir}/
├── Samples/{sample}/
│   ├── {sample}_unfiltered.h5ad                 CONTRACT (prepare_input_{vhd,x5k})
│   ├── QuPath_image/…                            generate_qupath_{vhd,x5k}
│   ├── validation/input_validation.json         validate_input
│   ├── qc_sweep/                                 qc_sweep (optional)
│   ├── adata_{sample}.h5ad, metadata_{sample}.tsv, {sample}_report.tsv   preprocess_umap
│   ├── leiden_resolution_*/                      leiden_analysis
│   ├── annotation/, adata_{sample}_annotated.h5ad   annotate_cells
│   └── neighbourhood_analysis/{annot_type}/      neighbourhood_analysis
├── integrated_samples/                           concatenated / harmony / sketched / samples_report.pdf
├── pseudobulk/{annot_type}/{analysis_level}/     aggregated/ + de_results/
├── Subcompartments/{name}/{Harmony,NoHarmony}/   subcluster
├── spatial_niches/                               tsv/, plots/
├── gene_exploration/                             expression_ranges.tsv, {entry}/…
└── cluster_annotations_template.tsv              generate_annotation_template
```

## Key design decisions

- **HeadBlock/CoreBlock split.** Technology-specific code produces a standardized unfiltered
  contract; the common CoreBlock consumes it. New platform = new HeadBlock, CoreBlock unchanged.
- **Coherent naming.** Head rules carry a 3-letter technology code (`_vhd`, `_x5k`) so the
  organization is easy to follow and heads can coexist.
- **`validate_input` is a DAG gate.** A hard failure raises → no report → the core cannot start
  for that sample. Soft issues (missing region annotation, no mito genes, no image) are recorded.
- **`qc_sweep` is diagnostic only** — it never filters, clusters, or writes an h5ad.
- **External annotation takes over.** When enabled it becomes the primary annotation everywhere,
  so cells the external tool left unannotated can't leak into downstream analyses (see
  [Configuration](configuration.md)).
- **Config-driven colours** — regions, sample metadata, and cell types, applied consistently across
  every plot, with a grey fallback.
- **Retries scale memory.** `mem_mb` grows with the attempt number, so an OOM-killed job is
  resubmitted with more RAM.

## Extending to a new platform

Write a new HeadBlock — a `rules/*_<code>.smk` chain plus a `<code>` entry in `KNOWN_HEADS` and a `mode` enum value — that emits the contract h5ad described above, then set `mode: <code>`. The entire CoreBlock runs unchanged. Alternatively, prepare contracts with any external tool and run the CoreBlock in `mode: decoupled`.
