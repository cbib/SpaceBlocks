This page explains how to fully set up the config files to run SpaceBlocks. Full documentation can be found at https://cbib.github.io/SpaceBlocks/configuration/. **If you are visualizing the documentation from the Snakemake catalog** we recommend visiting the SpaceBlocks GitHub, since the rendered tube map does not include the HeadBlocks and is not consistently colored.

Copy **`config/config.yaml.template`** to **`config/config.yaml`** before configuring a run. The template contains all parameters, while the committed `config/config.yaml` is only the runnable fixture used by the Snakemake Workflow Catalog.

The explanations on this page are divided by type (directory, parameter, color, etc.). In the config file, parameters are divided by Block and function.

**Every key in `config/config.yaml` is validated** against `workflow/schemas/config.schema.yaml` before the run starts, so a typo or a missing required field fails immediately with a clear message.

!!! warning "Complete list of SpaceBlocks parameters"
    The full, always-updated table of every parameter (type, default, required) is generated automatically from the schema and shown on the workflow's [Snakemake-catalog page](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/cbib/SpaceBlocks.html#workflow-parameters).
    This page covers the *how* and the *why*; the schema is the exhaustive reference.

## 1. Choose a mode

The config `mode` selects which **Headblock** builds the standardized contract h5ad.  In `decoupled` mode, SpaceBlocks consumes contract h5ads that already exist.

!!! warning "Space Ranger installation"
    The **Visium HD** HeadBlock additionally needs an external Space Ranger installation.

| `config["mode"]` | Head that runs | Input |
| --- | --- | --- |
| `visiumhd` | `spaceranger_count_vhd` → `generate_qupath_vhd` → `prepare_input_vhd` | 10x Visium HD fastqs + Space Ranger |
| `xenium5k` | `convert_zarr_x5k` → `generate_qupath_x5k` → `prepare_input_x5k` | Xenium output bundle |
| `atera` | `convert_zarr_ate` → `generate_qupath_ate` → `prepare_input_ate` | Atera (alpha) output bundle; optional registered H&E inputs |
| `merscope` | `generate_qupath_mer` → `prepare_input_mer` | MERSCOPE region directory |
| `decoupled` | *(none)* | Pre-existing contract h5ad files you provide |

The analysis **CoreBlocks** are identical in all five cases. Atera support is currently in alpha because the public preview format may change before commercial release.

## 2. Paths, files and sample sheets

The SpaceBlocks `config/config.yaml` sets the input and output directories, and gets the sample information and metadata (for integration and plotting) from one or two TSV files, depending on the `config["mode"]` set.

!!! tip "Setting a base directory"
    `base_dir` is the project-root prefix for generated output directories. In the shipped template, most output folders are written relative to it, but several required inputs are intentionally left as explicit filesystem paths because they point to external software or reference data.

Paths in SpaceBlocks fall into three categories:

- **Project output paths**: typically relative to `base_dir` or kept under the config directory for metadata files.
- **Project input paths**: configuration or pre-computed data used by the pipeline. These can be defined relative to `config/` or as absolute paths.
- **External absolute paths**: installed tools, references, and other files that must be supplied explicitly.

### Project output paths

These are the paths that are normally defined relative to `base_dir`.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| Path | `post_processing_outdir` | **Mandatory** | Output directory root; it will contain the result folders for the run. |
| Path | `logdir` | **Mandatory** | Directory for per-rule logs and benchmarks. |
| Path | `spaceranger_processing_outdir` | Head modes | Directory where the head writes heavy intermediates such as Space Ranger outputs or SpatialData Zarr stores. |

### Project input paths

These paths contain metadata, configuration or pre-processed results that are ingested by the pipeline. We provide some defaults under the `config/config.yaml.template`. Others are defined as placeholders `/path/to/dir/` and `/path/to/file.ext`, so make sure to configure them if needed.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| File | `core_samples` | **Mandatory** | Sample metadata sheet. The shipped template keeps this in `config/` as project metadata. |
| Path | `contract_dir` | `mode: decoupled` | Directory holding the pre-existing contract h5ads to analyse. |
| Path | `xenium5k.xenium_dir` | `mode: xenium5k` | Directory pattern for each Xenium bundle. |
| Path | `atera.atera_dir` | `mode: atera` | `{sample}` path pattern to each Atera `outs/` bundle. |
| Path | `merscope.merscope_dir` | `mode: merscope` | `{sample}` path pattern to each MERSCOPE region directory. |
| Path | `geojson_path` | **Optional, recommended** | Directory containing the GeoJSONs that annotate spatial regions for each sample (see the [QuPath tutorial](https://cbib.github.io/SpaceBlocks/qupath-tutorial/)).  In `decoupled` mode this directory is not searched; region labels must already be in the contract h5ad. |
| Path | `precomputed_metadata_dir` | **Optional** (reproducibility) | Directory of precomputed per-sample metadata TSVs. The metadata can contain precomputed clusters and/or external annotations. |
| Path | `spatial_niches.niche_dir` | **Optional** (reproducibility) | Directory for precomputed niche TSVs when `use_precomputed` is enabled. |
| File | `per_sample_qc` | **Optional** (per-sample filters) | TSV of sample-specific QC filters. If empty, the same `analysis.*` thresholds apply to every sample (*default*). |
| File | `snakemake_cell_markers` | **Optional** (pre-annotation) | TSV of canonical marker genes drawn on the Leiden diagnostic plots (`leiden_analysis`) to guide manual annotation. If not set, default markers are used (`config/snakemake_cell_markers.tsv`) |
| File | `gene_exploration.queries` | Exploration Coreblock | Genes / gene sets to score (AUCell) and plot in the exploration block. If not set, test gene queries are used (`config/gene_queries.tsv`). |
| File | `cluster_annotations` | **Optional** | TSV of cluster-to-cell-type labels used during annotation. Template ships with a template of placeholder annotations (`config/cluster_annotations.tsv`). If not set, cell will be labeled `Unannotated` |

!!! note "GeoJSON filenames"
    Region annotation is optional per sample, but a sample-prefixed GeoJSON with an unsupported suffix is treated as an error. Expected suffixes are `_tissue_hires_image.geojson` for Visium HD, `_morphology.geojson` for Xenium and MERSCOPE, and either `_he_background.geojson` or `_morphology.geojson` for Atera.

### Explicit external inputs

These are not derived from `base_dir` and must be set explicitly to the correct filesystem location.

| Path / File | Key | Condition | Meaning |
| --- | --- | --- | --- |
| File | `spaceranger` | `mode: visiumhd` | Path to the externally installed Space Ranger executable. |
| File | `probe_set` | `mode: visiumhd` | [Probe set CSV](https://www.10xgenomics.com/support/spatial-gene-expression-hd/documentation/steps/probe-sets) needed to run Space Ranger. |
| File | `transcriptome` | `mode: visiumhd` | [Reference transcriptome](https://www.10xgenomics.com/support/software/space-ranger/downloads) needed to run Space Ranger. |
| File | `ingest_ref` | **Optional** (auto-annotation) | Annotated h5ad scRNA-seq reference for `ingest_ref`. |
| File | `samples` | `mode: visiumhd` | Sample sheet with fastq directories, slide, and capture area per sample. |

### Sample sheets

`core_samples.tsv` is the technology-agnostic sample sheet used by the CoreBlocks in every use case. This file must contain a column named `sample`. Every other column is stamped into each sample's `obs`, so it can be used by integration plots and pseudobulk `group_by`, `covariates`, or `paired_by` settings

Visium HD additionally uses `visiumhd_samples.csv` for fastq, slide, and capture-area information.

For Xenium, Atera, and MERSCOPE, the sample list comes from `core_samples`, and each input bundle is located with the corresponding `{sample}` path pattern.

!!! tip "Color customization for result visualization"
    You can customize the color scale for any sample metadata in `core_samples.tsv` as additional columns. See [section 3](#3-parameters) and [section 4](#4-color-scale-customization) for details.

| Sample sheet | Condition | Role |
| --- | --- | --- |
| `config/core_samples.tsv` | **Mandatory** | The technology-agnostic anchor: one row per sample. The first column is the sample name; any further columns are design metadata (e.g. `patient`, `condition`, `type`) that are stamped into `obs` and can be surfaced downstream. |
| `config/visiumhd_samples.csv` | `mode: visiumhd` | The Visium HD head sheet: fastq directories (plus any re-sequencing runs), slide, and area per sample. |

## 3. Common parameters

Parameters are settings that live exclusively in `config/config.yaml` and determine the details of the run.

| Parameter | Condition | Role |
| --- | --- | --- |
| `mode` | **Mandatory** | Which headblock builds the contract (see [section 1](#1-choose-a-mode)). |
| `random_seed` | Reproducibility | Seed for stochastic steps, including PCA, UMAP, Leiden, subclustering, neighbours, and sketching. Results may still differ across software stacks or systems (see [section 6](#6-advanced-reproducibility-and-reusability)). |
| `use_precomputed_clusters` | **Optional** (reproducibility) | If `true`, reuse Leiden clusters/metadata from `precomputed_metadata_dir` instead of recomputing them. |
| `ingest_ref_label_key` | Annotated reference | Column in the `ingest_ref` reference that holds the reference cell-type labels. |
| `integration.integrate_key` | **Mandatory** | Variable Harmony corrects over during integration (e.g. `sample`). |
| `annotation_types` | `[external_annotation]` when external annotation is enabled; otherwise `[tsv_annotation]` | Annotation variants expanded by pseudobulk and neighbourhood analyses. Allowed values are `tsv_annotation`, `external_annotation`, and `ingest_annotation`. Add `ingest_annotation` explicitly if those analyses should use transferred labels. |
| `extra_annotations.columns` | **Optional** | `core_samples.tsv` columns to be used in downstream plots (e.g. `[patient, batch]`). All sample-sheet columns are carried into `obs` and can be used by pseudobulk models whether or not they are plotted. |
| `analysis.run_leiden_analysis` | **Optional**, default to `true` | Controls the per-resolution diagnostic plot rule. Leiden clusters are still computed by `preprocess_umap`. |
| `analysis.run_pseudobulk_de` | **Optional** `false` in runtime/schema; `true` in the template | Enables DESeq2 differential expression. Pseudobulk aggregation may still run when DE is disabled. |
| `analysis.pseudobulk.analyses` | Pseudobulk | Named aggregation/model specifications containing the comparison variables, explicit contrasts, exclusions, optional pairing/covariates, one-vs-rest biomarkers, and optional LRT. |
| `analysis.umap_point_size` | **Optional** | Marker size for ordinary UMAP embedding plots (default: `2`). Highlight and split UMAPs remain at least size `10`. |
| `analysis.spatial_point_size` | **Optional** | Marker size for spatial maps (default: `20`). This is independent of the UMAP marker size. Coordinate-only spatial-niche maps are scaled proportionally. |

## 4. Color scale customization

In SpaceBlocks, color scales are fully customizable and consistent across analysis plots.

!!! warning "Undefined variable levels are colored in grey"
    When customizing color palettes, levels not explicitly listed fall back to grey `#cccccc`, so a value that renders grey usually means a missing key. If you wish to customize visualization, list all levels in a variable and assign a color for each of them.

To set a given color scale, you just need to specify the HEX color code for each level in `config/config.yaml`. There are three palette families, each keyed by *column → level → hex*:

```yaml
# core_samples.tsv columns: column -> level -> colour
sample_colors:
  patient:
    "Patient 1": "#000000"
    "Patient 2": "#E69F00"
  condition:
    Control: "#0072B2"
    Treatment: "#D55E00"
  batch:
    "Batch 1": "#8E44AD"

# annotation obs columns: column -> level -> colour
annotation_colors:
  cell_type_tsv:
    Tcells: "#1f77b4"
    Fibroblasts: "#2ca02c"
  cell_type_external:       # not inherited from cell_type_tsv
    Tcells: "#1f77b4"
  cell_type_ingest:
    Tcells: "#1f77b4"
  spatial_niche:            # optional override of the deterministic niche palette
    "0": "#4C78A8"

analysis:
  # Preferred display order. For a single region group in pseudobulk,
  # this is also the fallback when condition_order is absent.
  region_levels: ["Tumor area", "Healthy area"]
  # Region level -> colour
  region_colors:
    "Tumor area": "#46337EFF"
    "Healthy area": "#FDE725FF"
```

Spatial niches are coloured automatically from a deterministic palette unless you add a `spatial_niche` block under `annotation_colors`.

`region_levels` controls ordering, but it does not filter observations. For pseudobulk analyses, you may use `exclude_levels` to remove values.

## 5. Key configuration sections

The rest of the configuration lives in comment-separated config chunks. Files and single parameters are covered in [section 2](#2-paths-files-and-sample-sheets) and [section 3](#3-parameters); colour palettes in [section 4](#4-color-scale-customization).

| Section | Condition | What it configures |
| --- | --- | --- |
| `visiumhd` | `mode: visiumhd` | Visium HD head settings: `spaceranger`, `probe_set`, and `transcriptome`. |
| `xenium5k` | `mode: xenium5k` | Xenium head settings: `xenium_dir`, pyramid levels, and `pixel_size_um`. |
| `atera` | `mode: atera` | Atera head settings: `atera_dir`, Zarr and pyramid options, plus optional registered H&E image, alignment, and keypoint patterns. |
| `merscope` | `mode: merscope` | MERSCOPE head settings: `merscope_dir`, selected z-plane, embedded-image resolution, and image channels. |
| `contract` | **Optional**, defaults are implemented | Contract keys such as `sample_key`, `spatial_key`, `require_region`, `require_raw_counts`, and `mito_prefix`. |
| `qc_sweep` | **Optional** | Candidate-threshold QC diagnostics (never filters). |
| `external_annotation` | **Optional** | Overlay labels from an external tool: `enabled`, `column`, `keep_unannotated` (see [section 6](#6-advanced-reproducibility-and-reusability)). |
| `analysis` | **Required** | The bulk of the run: QC filters (`min_counts` / `min_genes` / `min_cells` / `max_counts` / `max_pct_mt`), the Leiden `resolution_scan_*`, named pseudobulk analyses + thresholds, and the `run_*` toggles. Per-sample QC overrides come from `per_sample_qc`. |
| `spatial_niches` | **Optional** | BANKSY niche detection across concatenated samples. |
| `subcompartments` | **Optional** | Named cell-type subsets to re-cluster in `subcluster`. |
| `gene_exploration` | Exploration CoreBlock | Genes / gene sets to score (AUCell) and plot in the exploration block, plus its `niche_column` and rank fraction. |
| `annotation_colors` | Color customization | Annotations for samples, conditions, regions and other variables. |
| `resources` | **Required** | Per-rule `mem_mb` / `runtime` / `threads` (with a `default`). Memory scales with the retry attempt, so an OOM-killed job is resubmitted with more RAM. |

!!! note
    A few one-key blocks are documented elsewhere for readability: `integration` (`integrate_key`) and `extra_annotations` (`columns`) are parameters in [section 3](#3-parameters); `cluster_annotations` is a file in [section 2](#2-paths-files-and-sample-sheets); and the colour blocks (`sample_colors`, `annotation_colors`, `analysis.region_colors`) are in [section 4](#4-color-scale-customization).

### Spatial niches

`spatial_niches.enabled: true` runs BANKSY jointly across the preprocessed samples. The physical-neighbour and clustering-neighbour parameters are independent:

| Key | Default | Meaning |
| --- | ---: | --- |
| `lambda` | `0.8` | Weight of spatial-neighbour features in the BANKSY matrix. |
| `num_neighbours` | `18` | Physical neighbours used to construct BANKSY features. |
| `max_m` | `1` | Highest spatial-neighbourhood harmonic used by BANKSY. |
| `nbr_weight_decay` | `scaled_gaussian` | BANKSY neighbour-weight decay. |
| `use_hvg` | `false` | Restrict BANKSY to highly variable genes. The shipped configs set this to `true`. |
| `n_top_genes` | `analysis.n_top_genes`, otherwise `2000` | **Strongly recommended** due to heavy RAM usage. Number of HVGs when `use_hvg` is enabled. |
| `cluster_n_neighbours` | `50` | Neighbours in the separate Harmony-corrected BANKSY PCA graph used for clustering. |
| `niche_resolution` | `0.1` | Leiden resolution controlling niche granularity. |
| `n_iterations` | `-1` | Leiden iterations; `-1` runs until convergence and may take longer on very large datasets. |
| `compute_umap` | `false` | Whether to compute the optional integrated BANKSY UMAP. |
| `refine` | `false` | Apply spatial majority-vote refinement to computed or reused labels. |

The BANKSY matrix is dense and has approximately `(2 + max_m) × n_genes` features. For large datasets, use `use_hvg: true` and tune `n_top_genes` to control memory.

With `use_precomputed: true`, SpaceBlocks skips BANKSY/Leiden and reloads a complete `niche_<sample>.tsv` for every sample. Missing files or cells will prevent the rule from running.

### Pseudobulk experimental designs

Pseudobulk analyses are configured under `analysis.pseudobulk.analyses`. Each analysis has a stable `name`, one `aggregation`, one or more `group_by` columns, optional `covariates` and `paired_by`, exclusions, explicit pairwise contrasts, and optional one-vs-rest and LRT tests.

| `aggregation` | One pseudobulk observation per | Separate result set per |
| --- | --- | --- |
| `all_cells` | sample | — |
| `by_celltype` | sample and cell type | cell type |
| `by_region` | sample and region | — |
| `by_celltype_region` | sample, cell type and region | cell type |

`cell_type` is a portable alias: SpaceBlocks resolves it to the cell-type column for
the active `annotation_types` entry, so manual TSV, external, and ingest-transferred
cell-type annotations can use the same pseudobulk configuration.

`condition_order` optionally sets the displayed order of comparison-group labels. With multiple `group_by` columns, use the combined labels written to pseudobulk metadata, for example `condition=Control | treatment=Drug`. Unlisted observed levels are appended.

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
The exclusions are also applied to one-vs-rest tests and the LRT.

Here P1–P6 are independent samples, so `paired_by` is omitted. Use covariates only when the design supports estimating them without confounding. Use `aggregation: all_cells` for one whole-sample result rather than a result per cell type.

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

This tests one combination against another. It is **not** a formal interaction or difference-of-differences test.

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
  one_vs_rest:
    enabled: true
  lrt:
    enabled: true
```

For a paired Wald contrast, SpaceBlocks retains only matched units containing both
requested levels and interprets `min_replicates` as the minimum number of complete pairs.

For separate samples matched by patient, add the shared patient identifier to `core_samples.tsv` and set `paired_by: patient`, as in this example:

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

#### Repeated observations without pairing

It is generally a good idea to used paired contrasts when observations can be considered repeated measurements from the same experimental unit. Thus, if no `paired_by` is configured and a sample contributes pseudobulks to more than one comparison level, DE for that pooled or cell-type subgroup is aborted by default because those observations are not independent. The script records the reason in that subgroup's `ERROR.txt`, so it does not silently omit the repeated samples.

For an explicitly exploratory analysis, set `allow_repeated_unpaired: true` and leave `paired_by` unset. The model then treats the repeated pseudobulks as independent, ignores within-sample correlation, and may underestimate standard errors. When repeated observations are detected, the analysis writes `WARNING_repeated_unpaired.txt`. `paired_by` and `allow_repeated_unpaired: true` cannot be combined.

#### Multiple level designs

Multiple level designs unlock the possibility of (1) running likelihood ratio tests (LRT) for holistic comparisons across experimental levels and (2) one-vs-rest Wald tests to identify potential gene markers for each of the levels. SpaceBlocks provides native implementation of these statistical tests, explained below.

When `lrt.enabled` is true and at least three eligible comparison-group levels exist, SpaceBlocks runs an omnibus LRT test (using DESeq2 R implementation) followed by gene clustering using the DEGpatterns R package. This is carried out only when sufficient level and observations remain after exclusions, replicate filtering, and pairing. This is performed following the tutorial from [Harvard Chan Bioinformatics Core workshop](https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html).
With `paired_by`, the full model includes the paired unit term and the reduced model removes
only the comparison group. With fewer than three eligible levels, LRT and DEGpatterns are
skipped for that subgroup or cell type, and the reason is written to the pseudobulk rule log and
`LRT/SKIPPED_insufficient.txt`. When the LRT does run and finds enough significant genes, DEGpatterns groups their expression profiles. Note that no order or trend direction is imposed by the configuration. Results are written under `pseudobulk/.../LRT/`.

When `one_vs_rest.enabled` is true, SpaceBlocks also tests each eligible level against
all other eligible levels combined. It uses the same exclusions explained for the LRT test.
Results are written under `pseudobulk/.../one_vs_rest/`, including a TSV, volcano plot, and heatmaps for
each level. `unique_markers.tsv` contains genes significantly upregulated for one level but
not for any other one-vs-rest comparison; the accompanying barplot and heatmap summarize
these level-specific biomarkers. This is independent of the omnibus LRT and configured
pairwise contrasts.

## 6. Advanced: Reproducibility and reusability

Even though we have ensured the highest reproducibility standards when creating SpaceBlocks, some steps are never 100% reproducible between systems (e.g. UMAP calculation, Leiden clustering).

We therefore provide features for minimal file sharing/storage that tighten the reproducibility gap.

SpaceBlocks allows you to input:

- **Externally assembled h5ad AnnData objects** — run `mode: decoupled` and point `contract_dir` at the directory of pre-built contract h5ads. The heads are skipped; the core validates and analyses them directly.
- **Pre-computed clusters / annotations** — set `use_precomputed_clusters: true`  and provide `metadata_<sample>.tsv` files in `precomputed_metadata_dir`.
- **Precomputed niches:** set `spatial_niches.enabled: true`, `spatial_niches.use_precomputed: true`, and optionally `spatial_niches.niche_dir`.
- **External cell labels:** enable `external_annotation` as described below.

### External annotation

```yaml
external_annotation:
  enabled: true
  column: manual_celltype
  keep_unannotated: true
```

For every sample, provide `column` in `precomputed_metadata_dir/metadata_<sample>.tsv`. **In `decoupled` mode only**, the source column may instead be stored directly in each contract h5ad's `obs` when no external metadata directory is configured. SpaceBlocks normalizes the source to `cell_type_external` downstream.

- With `keep_unannotated: true`, normal QC is applied and unmatched cells remain as `Unannotated`. A mixed real/`Unannotated` column remains active; a placeholder-only column is omitted from annotation plots.
- With `keep_unannotated: false`, preprocessing retains only cells with real external labels and skips the configured cell-QC thresholds. Genes expressed in no retained cell are still removed.

When external annotation is enabled, it becomes the primary annotation used by reports, integration summaries, and gene exploration. If `ingest_ref` is also configured, both `cell_type_external` and `cell_type_ingest` are retained. Pseudobulk and neighbourhood analyses use only the variants listed in `annotation_types`; they do not automatically add ingest labels.

`annotate_cells` can operate without a manual cluster mapping and will set `cell_type_tsv` to `Unannotated`. However, the current default-target assembly still uses a readable, non-empty `cluster_annotations` file as the gate for annotation-dependent downstream targets. Keep such a file configured when running `rule all`, even for an external-primary analysis.


## 7. Mode-specific examples

The commented `config/config.yaml.template` is the full template. Here, these snippets show only the keys that distinguish each mode, and you may combine them with the common settings from the template.

**Visium HD**

```yaml
mode: visiumhd
samples: config/visiumhd_samples.csv
spaceranger: /path/to/spaceranger
probe_set: /path/to/probe_set.csv
transcriptome: /path/to/refdata-gex
spaceranger_processing_outdir: "{base_dir}/spaceranger_results"
```

**Xenium 5K**

```yaml
mode: xenium5k
xenium5k:
  xenium_dir: "/path/to/xenium/{sample}"
  qupath_pyramid_level: 3
  hires_pyramid_level: 3
  pixel_size_um: 0.2125
```

**Atera (alpha support)**

```yaml
mode: atera
atera:
  atera_dir: "/path/to/atera/{sample}/outs"
  qupath_pyramid_level: 3
  hires_pyramid_level: 3
  pixel_size_um: 0.2125
  he_image: ""       # optional {sample} pattern; set with he_alignment
  he_alignment: ""   # optional {sample} pattern; set with he_image
  he_keypoints: ""   # optional {sample} pattern
```

**MERSCOPE**

```yaml
mode: merscope
merscope:
  merscope_dir: "/path/to/merscope/{sample}"
  z_index: 3
  hires_pixel_size_um: 1.0
  channels: [DAPI, PolyT, Cellbound1, Cellbound2, Cellbound3]
```

**Pre-built contracts**

```yaml
mode: decoupled
contract_dir: /path/to/contract_h5ads     # one <sample>.h5ad file per sample
```

## 8. Minimal decoupled example

The following is a schema-valid structural skeleton. Replace the paths and resource values for your environment; each sample listed in `core_samples` must have `<contract_dir>/<sample>.h5ad`.

```yaml
mode: "decoupled"
core_samples: "config/core_samples.tsv"
contract_dir: "/path/to/contract_h5ads"
geojson_path: "/path/to/geojson"
precomputed_metadata_dir: "/path/to/precomputed_metadata/"
post_processing_outdir: "{base_dir}/results"
# ... see config/config.yaml.template for the full, commented template.
```

For a normal run, start by copying `config/config.yaml.template` to `config/config.yaml`, then adjust the paths and sample sheets to your data.

You may next read the [get started](https://cbib.github.io/SpaceBlocks/getting-started/) documentation and the [public data end-to-end example runs](https://cbib.github.io/SpaceBlocks/demos/).
