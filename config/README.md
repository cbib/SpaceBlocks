This page explains how to fully set up the config files to run SpaceBlocks. Full documentation can be found at https://cbib.github.io/SpaceBlocks/configuration/.

SpaceBlocks is configured through **`config/config.yaml`**, which contains all the parameters needed for the run, plus one or two **sample sheets**.

The explanations on this page are divided by type (directory, parameter, color, and so on). In the config file, parameters are divided by Block and function.

**Every key in `config/config.yaml` is validated** against `workflow/schemas/config.schema.yaml` before the run starts, so a typo or a missing required field fails immediately with a clear message.

!!! note "Snakemake-catalog page"
    The full, always-updated table of every parameter (type, default, required) is generated automatically from the schema and shown on the workflow's [Snakemake-catalog page](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/cbib/SpaceBlocks.html#workflow-parameters).
    This page covers the *how* and the *why*; the schema is the exhaustive reference.

## 1. Choose a mode

`mode` selects which **Headblock** builds the standardized contract h5ad (or none):

!!! warning "Space Ranger installation"
    The **Visium HD** HeadBlock additionally needs an external Space Ranger installation.

| `config["mode"]` | Head that runs | Input |
| --- | --- | --- |
| `visiumhd` | `spaceranger_count_vhd` → `generate_qupath_vhd` → `prepare_input_vhd` | 10x Visium HD fastqs + Space Ranger |
| `xenium5k` | `convert_zarr_x5k` → `generate_qupath_x5k` → `prepare_input_x5k` | Xenium output bundle |
| `atera` | `convert_zarr_ate` → `generate_qupath_ate` → `prepare_input_ate` | Atera output bundle; optional registered H&E inputs |
| `merscope` | `generate_qupath_mer` → `prepare_input_mer` | MERSCOPE region directory |
| `decoupled` | *(none)* | Pre-existing contract h5ad files you provide |

The analysis **CoreBlocks** are identical in all five cases. Atera support is currently alpha because the public preview format may change before commercial release.

## 2. Paths, files and sample sheets

The SpaceBlocks `config/config.yaml` sets the input and output directories, and gets the sample information and metadata (for integration and plotting) from one or two TSV files, depending on the `config["mode"]` set.

!!! tip "Setting a base directory"
    `base_dir` is the project-root prefix for generated output directories. In the shipped template, most output folders are written relative to it, but several required inputs are intentionally left as explicit filesystem paths because they point to external software or reference data.

Paths in SpaceBlocks fall into three categories:

- Project output paths: typically relative to `base_dir` or kept under the config directory for metadata files.
- Project input paths: configuration or pre-computed data used by the pipeline. These can be defined relative to `config/` or as absolute paths.
- External absolute paths: installed tools, references, and other files that must be supplied explicitly.

### Project output paths

These are the paths that are normally defined relative to `base_dir`.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| Path | `post_processing_outdir` | **Mandatory** | Output directory root; it will contain the result folders for the run. |
| Path | `logdir` | **Mandatory** | Directory for per-rule logs and benchmarks. |
| Path | `spaceranger_processing_outdir` | Head modes | Directory where the head writes heavy intermediates such as Space Ranger outputs or SpatialData Zarr stores. |

### Project input paths

These paths contain metadata, configuration or pre-processed results that are ingested by the pipeline. We provide some defaults under the `config/` directory. Others are defined as placeholders `/path/to/dir/` and `/path/to/file.ext`, so make sure to configure them if needed.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| File | `core_samples` | **Mandatory** | Sample metadata sheet. The shipped template keeps this in `config/` as project metadata. |
| Path | `contract_dir` | `mode: decoupled` | Directory holding the pre-existing contract h5ads to analyse. |
| Path | `xenium5k.xenium_dir` | `mode: xenium5k` | Directory pattern for each Xenium bundle. |
| Path | `spatial_niches.niche_dir` | `spatial_niches` | Directory for precomputed niche TSVs when `use_precomputed` is enabled. |
| Path | `geojson_path` | **Optional, recommended** | Directory containing the GeoJSONs that annotate spatial regions for each sample (see the [QuPath tutorial](https://cbib.github.io/SpaceBlocks/qupath-tutorial/)). |
| Path | `precomputed_metadata_dir` | **Optional** (reproducibility) | Directory of precomputed per-sample metadata TSVs. The metadata can contain precomputed clusters and/or external annotations. |
| File | `per_sample_qc` | **Optional** (per-sample filters) | TSV of sample-specific QC filters. If empty, the same `analysis.*` thresholds apply to every sample (*default*). |
| File | `snakemake_cell_markers` | **Optional** (pre-annotation) | TSV of canonical marker genes drawn on the Leiden diagnostic plots (`leiden_analysis`) to guide manual annotation. If not set, default markers are used (`config/snakemake_cell_markers.tsv`) |
| File | `gene_exploration.queries` | Exploration Coreblock | Genes / gene sets to score (AUCell) and plot in the exploration block. If not set, test gene queries are used (`config/gene_queries.tsv`). |
| File | `cluster_annotations` | **Optional** | TSV of cluster-to-cell-type labels used during annotation. Template ships with a template of placeholder annotations (`config/cluster_annotations.tsv`). If not set, cell will be labeled `Unannotated` |

### External absolute paths

These are not derived from `base_dir` and must be set explicitly to the correct filesystem location.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| File | `spaceranger` | `mode: visiumhd` | Path to the external Space Ranger installation. |
| File | `probe_set` | `mode: visiumhd` | [Probe set CSV](https://www.10xgenomics.com/support/spatial-gene-expression-hd/documentation/steps/probe-sets) needed to run Space Ranger. |
| File | `transcriptome` | `mode: visiumhd` | [Reference transcriptome](https://www.10xgenomics.com/support/software/space-ranger/downloads) needed to run Space Ranger. |
| File | `ingest_ref` | **Optional** (auto-annotation) | Annotated scRNA-seq reference h5ad for `ingest_ref`. |
| File | `samples` | `mode: visiumhd` | Sample sheet with fastq directories, slide, and capture area per sample. |

The template uses placeholder values such as `/path/to/...` for these external inputs to make it clear they must be replaced before running the workflow.

`core_samples.tsv` is the technology-agnostic sample sheet used by the CoreBlocks in every use case. Visium HD additionally uses `visiumhd_samples.csv` for fastq, slide, and capture-area information.

For `xenium5k`, `atera`, and `merscope`, per-sample input bundles are located through the corresponding `{sample}` path pattern in `xenium5k.xenium_dir`, `atera.atera_dir`, or `merscope.merscope_dir`; no second head-specific sample sheet is required.

!!! tip "Color customization for result visualization"
    You can customize the color scale for any sample metadata in `core_samples.tsv` as additional columns. See [section 3](#3-parameters) and [section 4](#4-color-scale-customization) for details.

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
| `extra_annotations.columns` | **Optional** | `core_samples.tsv` columns to surface in downstream plots (e.g. `[patient, batch]`). All sample-sheet columns are carried into `obs` and can be used by pseudobulk models whether or not they are plotted. |
| `analysis.pseudobulk.analyses` | Pseudobulk | Named aggregation/model specifications containing the comparison variables, explicit contrasts, exclusions, optional pairing/covariates, and optional LRT. |

## 4. Color scale customization

In SpaceBlocks, color scales are fully customizable and consistent across analysis plots.

!!! warning "Undefined variable levels are colored in grey"
    When customizing color palettes, levels not explicitly listed fall back to grey, so a value that renders grey usually means a missing key. If you wish to customize visualization, list all levels in a variable and assign a color for each of them.

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

# Preferred region order in QC-sweep, subcluster, and region-only pseudobulk split heatmaps.
# Entries are shown first-to-last (typically left-to-right); unlisted levels follow.
# region_colors sets the palette. Pseudobulk filtering uses exclude_levels.
analysis:
  region_levels: ["Tumor area", "Healthy area"]
  region_colors:
    "Tumor area": "#46337EFF"
    "Healthy area": "#FDE725FF"
```

Spatial niches are coloured automatically from a deterministic palette unless you add a `spatial_niche` block under `annotation_colors`.

## 5. Key sections

The rest of the configuration lives in nested blocks. Files and single parameters are covered in [section 2](#2-paths-files-and-sample-sheets) and [section 3](#3-parameters); colour palettes in [section 4](#4-color-scale-customization).

| Section | Condition | What it configures |
| --- | --- | --- |
| `xenium5k` | `mode: xenium5k` | Xenium head settings: `xenium_dir`, `zarr_dir`, pyramid levels, and `pixel_size_um`. |
| `atera` | `mode: atera` | Atera head settings: `atera_dir`, Zarr and pyramid options, plus optional registered H&E image, alignment, and keypoint patterns. |
| `merscope` | `mode: merscope` | MERSCOPE head settings: `merscope_dir`, selected z-plane, embedded-image resolution, and image channels. |
| `contract` | **Mandatory** | Semantic keys of the hand-off object: `sample_key`, `spatial_key`, `require_region`, `require_raw_counts`, `mito_prefix`. |
| `analysis` | **Mandatory** | The bulk of the run: QC filters (`min_counts` / `min_genes` / `min_cells` / `max_counts` / `max_pct_mt`), the Leiden `resolution_scan_*`, named pseudobulk analyses + thresholds, and the `run_*` toggles. Per-sample QC overrides come from `per_sample_qc`. |
| `resources` | **Mandatory** | Per-rule `mem_mb` / `runtime` / `threads` (with a `default`). Memory scales with the retry attempt, so an OOM-killed job is resubmitted with more RAM. |
| `external_annotation` | *optional* | Overlay labels from an external tool: `enabled`, `column`, `keep_unannotated` (see [section 6](#6-advanced-reproducibility-and-reusability)). |
| `qc_sweep` | *optional* | Candidate-threshold QC diagnostics (never filters). |
| `spatial_niches` | *optional* | BANKSY niche detection across concatenated samples. |
| `subcompartments` | *optional* | Named cell-type subsets to re-cluster in `subcluster`. |
| `gene_exploration` | Exploration CoreBlock | Genes / gene sets to score (AUCell) and plot in the exploration block, plus its `niche_column` and rank fraction. |

!!! note
    A few one-key blocks are documented elsewhere for readability: `integration` (`integrate_key`) and `extra_annotations` (`columns`) are parameters in [section 3](#3-parameters); `cluster_annotations` is a file in [section 2](#2-paths-files-and-sample-sheets); and the colour blocks (`sample_colors`, `annotation_colors`, `analysis.region_colors`) are in [section 4](#4-color-scale-customization).

### Pseudobulk experimental designs

Pseudobulk analyses are configured under `analysis.pseudobulk.analyses`. Each entry
has a stable `name`, one aggregation strategy, the sample metadata that define the
comparison groups, and the statistical tests to run. Counts are summed within each
pseudobulk observation described in the table below.

The available aggregations are:

| `aggregation` | One pseudobulk observation per | Separate result set per |
| --- | --- | --- |
| `all_cells` | sample | — |
| `by_celltype` | sample and cell type | cell type |
| `by_region` | sample and region | — |
| `by_celltype_region` | sample, cell type and region | cell type |

`cell_type` is a portable alias: SpaceBlocks resolves it to the cell-type column for
the active `annotation_types` entry, so manual TSV, external, ingest-transferred, and
refined cell-type annotations can use the same pseudobulk configuration.

#### Independent phenotype comparison

Add the experimental variables to `core_samples.tsv`:

```tsv
sample	phenotype	batch
Cancer_P1	Cancer	Batch1
Cancer_P2	Cancer	Batch2
Cancer_P3	Cancer	Batch3
Normal_P4	Normal	Batch1
Normal_P5	Normal	Batch2
Normal_P6	Normal	Batch3
```

Then compare phenotypes for each cell type while adjusting for batch:

```yaml
analysis:
  run_pseudobulk_de: true
  min_replicates: 3
  pseudobulk:
    analyses:
      - name: phenotype_by_celltype
        aggregation: by_celltype
        group_by: [phenotype]
        covariates: [batch]
        exclude_levels:
          cell_type: [Unannotated]
        contrasts:
          - name: cancer_vs_normal
            numerator: {phenotype: Cancer}
            denominator: {phenotype: Normal}
        lrt:
          enabled: false
```

The numerator determines the positive log2-fold-change direction. Levels not named
in a pairwise contrast can remain in the fitted model and help dispersion estimation.
Use `exclude_levels` only for values that should be removed from the analysis entirely.
The exclusions are also applied to the LRT.

In this example, P1–P6 are six independent biological samples: no individual
contributes to both phenotypes, so `paired_by` must be omitted. The batch covariate is
valid because every batch contains one Cancer and one Normal sample; omit `covariates`
when there is no nuisance variable to adjust for. Use `aggregation: all_cells` instead
if one whole-sample result is wanted rather than a separate result for each cell type.

#### Phenotype and treatment

Multiple `group_by` columns form a combined categorical group. In the example below,
Cancer and Normal phenotype levels are combined with Drug and Vehicle treatment
levels. This supports clear comparisons between observed combinations without exposing
interaction formulas:

```yaml
- name: treatment_by_celltype
  aggregation: by_celltype
  group_by: [phenotype, treatment]
  exclude_levels:
    cell_type: [Unannotated]
  contrasts:
    - name: drug_vs_vehicle_in_cancer
      numerator: {phenotype: Cancer, treatment: Drug}
      denominator: {phenotype: Cancer, treatment: Vehicle}
```

This tests one combination against another; it is **not** a formal statistical
interaction or difference-of-differences test.

#### Paired regions or matched samples

Use `paired_by: sample` when the same sample contributes both levels of a contrast,
as in tumour and healthy regions from the same section:

```yaml
- name: regions_by_celltype
  aggregation: by_celltype_region
  group_by: [region_annotation]
  paired_by: sample
  exclude_levels:
    region_annotation: [Unlabeled, Bubble]
    cell_type: [Unannotated]
  contrasts:
    - name: tumor_vs_healthy
      numerator: {region_annotation: Tumor area}
      denominator: {region_annotation: Healthy area}
  lrt:
    enabled: true
```

For a paired Wald contrast, SpaceBlocks retains only matched units containing both
requested levels and interprets `min_replicates` as the minimum number of complete pairs.

When Cancer and Normal are separate samples matched from the same patients, add a
shared patient identifier to `core_samples.tsv` and pair by that identifier:

```tsv
sample	phenotype	patient
Cancer_P1	Cancer	P1
Normal_P1	Normal	P1
Cancer_P2	Cancer	P2
Normal_P2	Normal	P2
Cancer_P3	Cancer	P3
Normal_P3	Normal	P3
```

```yaml
- name: matched_phenotype_by_celltype
  aggregation: by_celltype
  group_by: [phenotype]
  paired_by: patient
  contrasts:
    - name: cancer_vs_normal
      numerator: {phenotype: Cancer}
      denominator: {phenotype: Normal}
```

Thus, use `paired_by: sample` for multiple region levels from the same spatial sample,
`paired_by: patient` for separate samples matched by patient, and omit `paired_by` for
independent Cancer-versus-Normal samples.

When `lrt.enabled` is true, SpaceBlocks runs an omnibus test across every non-excluded
level with sufficient replication. With `paired_by`, the full model includes the paired
unit term and the reduced model removes only the comparison group. DEGpatterns is run
when the LRT finds enough significant genes; no order or trend direction is imposed by
the configuration.

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
spaceranger_processing_outdir: "{base_dir}/sr_out"
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

**Atera** (`mode: atera`, alpha)
```yaml
mode: "atera"
atera:
  atera_dir: "/path/to/atera/{sample}/outs"
  zarr_dir: ""
  qupath_pyramid_level: 3
  hires_pyramid_level: 3
  pixel_size_um: 0.2125
  he_image: ""       # optional {sample} pattern; set with he_alignment
  he_alignment: ""   # optional {sample} pattern; set with he_image
  he_keypoints: ""   # optional {sample} pattern
```

**MERSCOPE** (`mode: merscope`)
```yaml
mode: "merscope"
merscope:
  merscope_dir: "/path/to/merscope/{sample}"
  z_index: 3
  hires_pixel_size_um: 1.0
  channels: [DAPI, PolyT, Cellbound1, Cellbound2, Cellbound3]
```

**Externally prepared data** (`mode: decoupled`)
```yaml
mode: "decoupled"
contract_dir: "/path/to/contract_h5ads"    # one <sample>.h5ad file per sample
```

All modes additionally set `core_samples`, `post_processing_outdir`, `logdir`, and the other required common keys in `config/config.yaml`; **`geojson_path` is optional but strongly recommended** and may be left as a project-local directory or an explicit path.

## 8. Minimal example

```yaml
mode: "visiumhd"
samples: "config/visiumhd_samples.csv"
core_samples: "config/core_samples.tsv"
geojson_path: "path/to/geojson"
precomputed_metadata_dir: "/path/to/metadata_visiumhd/"
post_processing_outdir: "{base_dir}/results"
# ... see config/config.yaml for the full, commented template.
```

Start from the commented `config/config.yaml` shipped with the workflow and adjust the paths and sample sheets to your data.

You may next read the [get started](https://cbib.github.io/SpaceBlocks/getting-started/) documentation and the [public data end-to-end example runs](https://cbib.github.io/SpaceBlocks/demos/).
