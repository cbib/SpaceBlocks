# Configuration

This page explains how to fully set up the config files to run SpaceBlocks. Full documentation can be found at https://cbib.github.io/SpaceBlocks/configuration/.

SpaceBlocks is configured through **`config/config.yaml`**, which contains all the parameters needed for the run, plus one or two **sample sheets**.

The explanations in this page are divided by type (dir, parameter, color, etc.). In the config file, parameters are divided by Block and function.

**Every key in `config/config.yaml` is validated** against `workflow/schemas/config.schema.yaml` before the run starts, so a typo or a missing required field fails immediately with a clear message.

> [!NOTE]
> The full, always-updated table of every parameter (type, default, required) is generated automatically from the schema and shown on the workflow's Snakemake-catalog page.
> This page covers the *how* and the *why*; the schema is the exhaustive reference.

**On this page:**

1. [Choose a mode](#1-choose-a-mode)
2. [Paths, files and sample sheets](#2-paths-files-and-sample-sheets)
3. [Parameters](#3-parameters)
4. [Color scale customization](#4-color-scale-customization)
5. [Key sections](#5-key-sections)
6. [Reproducibility & reusability](#6-advanced-reproducibility-and-reusability)
7. [Example use cases](#7-example-use-case-configurations)
8. [Minimal example](#8-minimal-example)

## 1. Choose a mode

`mode` selects which **Headblock** builds the standardized contract h5ad (or none):

> [!WARNING]
> The **Visium HD** HeadBlock additionally needs an external Space Ranger installation.

| `config["mode"]` | Head that runs | Input |
| --- | --- | --- |
| `visiumhd` | `spaceranger_count_vhd` → `generate_qupath_vhd` → `prepare_input_vhd` | 10x Visium HD fastqs + Space Ranger |
| `xenium5k` | `convert_zarr_x5k` → `generate_qupath_x5k` → `prepare_input_x5k` | Xenium output bundle |
| `decoupled` | *(none)* | pre-existing contract h5ads you provide |

The analysis **Coreblock** is identical in all three cases.

## 2. Paths, files and sample sheets

The SpaceBlocks `config/config.yaml` sets the input and output directories, and gets the sample information and metadata (for integration and plotting) from one or two TSV files, depending on the `config["mode"]` set.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| Path | `spaceranger_processing_outdir` | `mode: visiumhd` / `xenium5k` | Directory where the head writes its heavy intermediates (Space Ranger tree / SpatialData zarr). |
| File | `spaceranger` | `mode: visiumhd` | Path to the Space Ranger installation. |
| File | `probe_set` | `mode: visiumhd` | [Probe set CSV](https://www.10xgenomics.com/support/spatial-gene-expression-hd/documentation/steps/probe-sets) needed to run Space Ranger. |
| File | `transcriptome` | `mode: visiumhd` | [Reference transcriptome](https://www.10xgenomics.com/support/software/space-ranger/downloads) needed to run Space Ranger. |
| Path | `post_processing_outdir` | **Mandatory** | Output directory root; it will contain the result folders for the run. |
| Path | `logdir` | **Mandatory** | Directory for per-rule logs and benchmarks. |
| Path | `geojson_path` | **Optional, recommended** | Directory containing the GeoJSONs that annotate spatial regions for each sample (see the [QuPath tutorial](https://cbib.github.io/SpaceBlocks/qupath-tutorial/)). |
| Path | `contract_dir` | `mode: decoupled` | Directory holding the pre-existing contract h5ads to analyse. |
| File | `per_sample_qc` | **Optional** (per-sample filters) | TSV of sample-specific QC filters. If empty, the same `analysis.*` thresholds apply to every sample (*default*). |
| File | `snakemake_cell_markers` | **Optional** (pre-annotation) | TSV of canonical marker genes drawn on the Leiden diagnostic plots (`leiden_analysis`) to guide manual annotation. |
| File | `gene_exploration.queries` (`config/gene_queries.tsv`) | Exploration Coreblock | Genes / gene sets to score (AUCell) and plot in the exploration block. |
| File | `ingest_ref` | **Optional** (auto-annotation) | Annotated scRNA-seq reference h5ad for `ingest_ref`. |
| Path | `precomputed_metadata_dir` | **Optional** (reproducibility) | Directory of precomputed per-sample metadata TSVs. The metadata can contain precomputed clusters and/or external annotations. |

You only need one sample sheet (`core_samples.tsv`) for `mode: decoupled` and `mode: xenium5k`, and two (`core_samples.tsv` and `visiumhd_samples.csv`) for `mode: visiumhd`.

For `xenium5k`, the per-sample input bundles come from `xenium5k.xenium_dir` (filenames with a `{sample}` pattern), not a separate sheet.

> [!TIP]
> You can customize the color scale for any sample metadata in `core_samples.tsv` as additional columns. See [section 3](#3-parameters) and [section 4](#4-color-scale-customization) for details.

| Sample sheet | Condition | Role |
| --- | --- | --- |
| `config/core_samples.tsv` | **Mandatory** | The technology-agnostic anchor: one row per sample. The first column is the sample name; any further columns are design metadata (e.g. `patient`, `condition`, `type`) that are stamped into `obs` and can be surfaced downstream. |
| `config/visiumhd_samples.csv` | `mode: visiumhd` | The Visium HD head sheet: fastq directories (plus any re-sequencing runs), slide, and area per sample. |

## 3. Parameters

Parameters are settings that live exclusively in `config/config.yaml` and determine the details of the run.

| Parameter | Condition | Role |
| --- | --- | --- |
| `mode` | **Mandatory** | Which headblock builds the contract (see [section 1](#1-choose-a-mode)). |
| `random_seed` | Reproducibility | Seed for stochastic steps (subsampling, sketching). Note that UMAP/Leiden are not fully deterministic across systems (see [section 6](#6-advanced-reproducibility-and-reusability)). |
| `use_precomputed_clusters` | **Optional** (reproducibility) | If `true`, reuse Leiden clusters/metadata from `precomputed_metadata_dir` instead of recomputing them. |
| `ingest_ref_label_key` | `ingest_ref` reference set | Column in the `ingest_ref` reference that holds the reference cell-type labels. |
| `integration.integrate_key` | **Mandatory** | Variable Harmony corrects over during integration (e.g. `sample`). |
| `extra_annotations.columns` | **Optional** | `core_samples.tsv` columns to carry into `obs` and surface in downstream plots (e.g. `[patient, batch]`). |

## 4. Color scale customization

In SpaceBlocks, color scales are fully customizable and consistent across analysis plots.

> [!WARNING]
> When customizing color palettes, levels not explicitly listed fall back to grey, so a value that renders grey usually means a missing key.

To set a given color scale, you just need to specify the HEX color code for each level in `config/config.yaml`. There are three palette families, each keyed by *column → level → hex*:

```yaml
# Sample-metadata palettes — one block per column listed in extra_annotations.columns
sample_colors:
  patient:
    "Patient 1": "#000000"
    "Patient 2": "#E69F00"
  batch:
    "Batch 1": "#8E44AD"

# Cell-type palettes — one block per annotation column (tsv / external / …)
annotation_colors:
  cell_type_tsv:
    "Tcells": "#1f77b4"
    "Fibroblasts": "#2ca02c"
  cell_type_external:          # declare explicitly; not inherited from cell_type_tsv
    "Tcells": "#1f77b4"

# Region palette — lives under analysis, keyed by the region_levels
analysis:
  region_colors:
    "Tumor area": "#46337EFF"
    "Healthy area": "#FDE725FF"
```

Spatial niches are coloured automatically from a deterministic palette unless you add a `spatial_niche` block under `annotation_colors`.

## 5. Key sections

The rest of the configuration lives in nested blocks. Files and single parameters are covered in [section 2](#2-paths-files-and-sample-sheets) and [section 3](#3-parameters); colour palettes in [section 4](#4-color-scale-customization).

| Section | Condition | What it configures |
| --- | --- | --- |
| `xenium5k` | `mode: xenium5k` | Xenium head settings: `xenium_dir`, `zarr_dir`, pyramid levels, `pixel_size_um`. |
| `contract` | **Mandatory** | Semantic keys of the hand-off object: `sample_key`, `spatial_key`, `require_region`, `require_raw_counts`, `mito_prefix`. |
| `analysis` | **Mandatory** | The bulk of the run: QC filters (`min_counts` / `min_genes` / `min_cells` / `max_counts` / `max_pct_mt`), the Leiden `resolution_scan_*`, the pseudobulk `analysis_levels` + thresholds, the `region_levels`, and the `run_*` toggles. Per-sample QC overrides come from `per_sample_qc`. |
| `resources` | **Mandatory** | Per-rule `mem_mb` / `runtime` / `threads` (with a `default`). Memory scales with the retry attempt, so an OOM-killed job is resubmitted with more RAM. |
| `external_annotation` | *optional* | Overlay labels from an external tool: `enabled`, `column`, `keep_unannotated` (see [section 6](#6-advanced-reproducibility-and-reusability)). |
| `qc_sweep` | *optional* | Candidate-threshold QC diagnostics (never filters). |
| `spatial_niches` | *optional* | BANKSY niche detection across concatenated samples. |
| `subcompartments` | *optional* | Named cell-type subsets to re-cluster in `subcluster`. |
| `gene_exploration` | Exploration CoreBlock | Genes / gene sets to score (AUCell) and plot in the exploration block, plus its `niche_column` and rank fraction. |

> [!NOTE]
> A few one-key blocks are documented elsewhere for readability: `integration` (`integrate_key`) and `extra_annotations` (`columns`) are parameters in [section 3](#3-parameters); `cluster_annotations` is a file in [section 2](#2-paths-files-and-sample-sheets); and the colour blocks (`sample_colors`, `annotation_colors`, `analysis.region_colors`) are in [section 4](#4-color-scale-customization).

## 6. Advanced: Reproducibility and reusability

Even though we have ensured the highest reproducibility standards when creating SpaceBlocks, some steps are never 100% reproducible between systems (e.g. UMAP calculation, Leiden clustering).

We therefore provide features for minimal file sharing/storage that tighten the reproducibility gap.

SpaceBlocks allows you to input:

- **Externally assembled h5ad AnnData objects** — run `mode: decoupled` and point `contract_dir` at the directory of pre-built contract h5ads. The heads are skipped; the core validates and analyses them directly.
- **Pre-computed clusters / annotations** — set `use_precomputed_clusters: true` and `precomputed_metadata_dir` to reuse Leiden clusters and metadata; for niches, set `spatial_niches.use_precomputed: true` with `spatial_niches.niche_dir`.
- **Externally annotated data** — set `external_annotation.enabled: true` with `external_annotation.column`, and choose whether to keep or discard unannotated barcodes downstream via `external_annotation.keep_unannotated` (see [section 5](#5-key-sections)).

These files are generated during the run, and can be shared with minimum effort to reproduce downstream results from raw data.

## 7. Example use case configurations

The commented `config/config.yaml` is the full template. The mode-specific keys that differ are:

**Visium HD** (`mode: visiumhd`)
```yaml
mode: "visiumhd"
samples: "config/visiumhd_samples.csv"          # fastq dirs / slide / area per sample
spaceranger: "/path/to/spaceranger"
probe_set: "/path/to/probe_set.csv"
transcriptome: "/path/to/refdata-gex"
spaceranger_processing_outdir: "/path/to/sr_out"
```

**Xenium 5K** (`mode: xenium5k`)
```yaml
mode: "xenium5k"
xenium5k:
  xenium_dir: "/path/to/xenium/{sample}"   # {sample} pattern to each bundle
  zarr_dir: ""                             # "" → spaceranger_processing_outdir
  qupath_pyramid_level: 3
  hires_pyramid_level: 3
  pixel_size_um: 0.2125
```

**Externally prepared data** (`mode: decoupled`)
```yaml
mode: "decoupled"
contract_dir: "/path/to/contract_h5ads"    # one {sample}_unfiltered.h5ad per sample
```

All three additionally set `core_samples`, `post_processing_outdir`, `logdir`, and (recommended) `geojson_path`.

## 8. Minimal example

```yaml
mode: "visiumhd"
samples: "config/visiumhd_samples.csv"
core_samples: "config/core_samples.tsv"
geojson_path: "/path/to/geojson"
post_processing_outdir: "/path/to/results"
# ... see config/config.yaml for the full, commented template.
```

Start from the commented `config/config.yaml` shipped with the workflow and adjust the paths and sample sheets to your data.

You may next read the [get started](../docs/getting-started.md) documentation and the [examples with public data](../docs/reproduction.md).
