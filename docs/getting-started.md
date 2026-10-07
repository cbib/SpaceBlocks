# Get started: a full run

This page walks through a complete SpaceBlocks run and the recommendations to make it smooth. For the exhaustive key reference see [Configuration](configuration.md); for what each rule does, [Rule reference](rules.md).

If you want to see an example run on public data, see the [public data end-to-end example runs](demos.md) page — Visium HD and Xenium 5K worked examples can be found there.

!!! tip "Double check your Snakemake run"
    A run is defined almost entirely by two things: your **experimental design** (`core_samples.tsv`) and your **cell-type annotations**. Everything else has sensible defaults. Keep those two under version control and your run will be reproducible end-to-end.

## 1. Prerequisites

- [Snakemake](https://snakemake.readthedocs.io) ≥ 9.13 and Conda/Mamba.
- For the **Visium HD** head only: an external [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger) ≥ 4.0.1.
- [QuPath](https://qupath.github.io/) (desktop), recommended to include manual region annotations for the downstream analyses (see [QuPath annotation tutorial](qupath-tutorial.md)).
- Use `--profile profiles/default` for local execution, or adapt the provided `--profile profiles/slurm` example for a Slurm HPC (tested with version 24.11.5).

We recommend to always run Snakemake with `--sdm conda`, so each rule is executed within its dedicated environment.

## 2. Configure

Create your run configuration, then pick a `mode` and fill the sample sheet:

```bash
cp config/config.yaml.template config/config.yaml
```

1. **`config/core_samples.tsv`** — one row per sample; extra columns (e.g. `patient`, `condition`) become design metadata, and can be used for downstream analyses and visualization.
2. **`config/config.yaml`** — set `mode`, the input and output paths, and (for Visium HD) the head files. See [Configuration](configuration.md) for full details.
3. If you are running `mode : visiumhd`, you will need to fill `visiumhd_samples.csv` as well.

Before launching SpaceBlocks, we recommend testing your configuration using a dry-run:

```bash
snakemake -n --sdm conda      # dry-run: builds the DAG and validates the config
```
This reports invalid configuration values, builds the DAG, and shows which rule-specific environments would be used. Environment creation happens during the first real run.


## 3. Region annotation (*optional*, but **recommended**)

When using any of the HeadBlocks, you may generate the images that will be included in the h5ad AnnData objects:

```bash
snakemake qupath_images --sdm conda
```
You may then manually annotate spatial regions in QuPath, export GeoJSONs and make them available for SpaceBlocks via `geojson_path` in `config/config.yaml`.

Then follow the [QuPath annotation tutorial](qupath-tutorial.md), using the filename for your HeadBlock and saving each GeoJSON into `config["geojson_path"]`.

!!! note  "Alternatives for region annotation"
    We provide instructions for GeoJSON generation via QuPath, but alternative software (i.e. Napari) could be used to generate the region annotation files.

While this step is optional (without it, all observations will become annotated as `Unlabeled` region), it is extremely powerful to discriminate between morphologically different areas of the tissue across your samples. We strongly encourage users to make use of this SpaceBlocks feature.

### 3.1. Region annotation in decoupled mode

In `mode: decoupled` there is no HeadBlock, so `qupath_images` produces nothing and
`geojson_path` is not consumed. Region labels must already be stored in each contract as
`obs["region_annotation"]`.

If that column is absent, validation emits a soft warning under the `/results/Samples/<sample>/validation/` directory, preprocessing labels every cell `Unlabeled`, and region-dependent plots and analyses are skipped.
Set `contract.require_region: true` if missing regions should instead stop the workflow.

!!! warning "Coordinate scale divergence"
    Make sure the pixel size and coordinate system of the image you annotate match the contract's `obsm["spatial"]`. In head modes the `generate_qupath_*` rules export a correctly-scaled image; in `mode: decoupled` you are responsible for annotating on an image consistent with the coordinates you embedded.

## 4. Preprocess, then look before you filter

SpaceBlocks provides a rule to screen how different QC filters (under `qc_sweep:` -> `min_genes`, `min_counts`, `max_counts`, `max_pct_mt`) would affect each of your samples. QC candidates can be configured under `qc_sweep.thresholds`.

You can combine this screening with an external annotated reference (`ingest_ref` and `ingest_ref_label_key` in the config) to observe how many, where in the tissue and which (automatically annotated) cell type would be removed when applying a threshold. To run it, use:

```bash
snakemake qc_sweep_all      --sdm conda   # OPTIONAL: see where each candidate threshold cuts
```

**We recommend running `qc_sweep_all` first** to choose QC thresholds from your own data rather than the defaults.
You can next use the same thresholds for all samples or use a `per_sample_qc` TSV if samples need different cut-offs.

After defining your general (*default*) or sample-specific thresholds, you should specify the preprocessing parameters (under `analysis:` in the config).

!!! tip "Include literature markers during preprocessing"
    We **strongly recommend** to include some literature predefined markers to screen for cell types present in your samples (defined via `snakemake_cell_markers`).

Next, you may run the preprocessing CoreBlock:

```bash
snakemake run_preprocessing --sdm conda   # validate → QC → normalise → cluster
```

The output (see [Outputs](outputs.md) for what every file is, and [the output tree](design.md#output-structure) for where it lands) will show spatial plots and UMAPs containing clustering metrics (clustree, Silhouette scores) for each resolution in the range, predefined markers (*optional*), QC plots to diagnose potentially noisy clusters. Then, you may proceed to the annotation step.

When `spatial_niches.enabled: true`, this target also runs BANKSY jointly across all
preprocessed samples. Its dense matrix can require substantially more memory than the
per-sample steps. We recommend to use HVGs and set the `spatial_niches` resource block accordingly.

## 5. Annotate cell types

Annotation is the most limiting step in every high-throughput single cell pipeline. It conditions downstream analysis and, thus, it is a breakpoint between CoreBlocks within the SpaceBlocks workflow.

Through the development of SpaceBlocks, we have tested several automatic annotation tools in Visium HD data. The results in our own data had very limited accuracy and, thus, manual annotation is the default option.

### 5.1 Manual annotation

SpaceBlocks automatically generates a cluster-to-cell type annotation TSV template within preprocessing. In this TSV, rows are cluster numbers and columns are sample names and Leiden resolution (i.e. "Sample1_res0.6").

The resolution is automatically read by the `annotate_cells` rule, so in the TSV you need to update:
1. The equivalence between the cluster number and the cell type present.
2. The resolution desired to perform annotation in the column names.

Note that you may use different clustering resolutions, as long as they are found in the metadata, to annotate each of the samples. Change the header of the TSV template to your desired `{sample}_{resolution}` to do so.

!!! danger "Save a backup of your metadata files"
    Leiden clustering between two different runs may vary slightly, for example, changing the order of specific clusters. To avoid wasting time re-annotating your samples, we recommend users to back up the metadata files from the results to be annotated. They can be found under `{base_dir}/results/Samples/<sample>/metadata_<sample>.tsv`.

After filling the annotation TSV, point `config["cluster_annotations"]` at it. This will unlock postprocessing rules.

### 5.2 External annotation

You may use labels generated outside SpaceBlocks without a manual cluster-to-cell-type TSV.
Enable `external_annotation`, set its source `column`, and provide one
`metadata_{sample}.tsv` per sample under `precomputed_metadata_dir`. These files must be
indexed by the same cell identifiers as the contract. In `mode: decoupled`, the configured
column may instead be stored directly in each contract's `obs`.

With `keep_unannotated: true`, ordinary QC is applied and unmatched cells remain as
`Unannotated`. With `keep_unannotated: false`, preprocessing keeps only cells with a real
external label and skips the ordinary QC thresholds because the external annotation defines
the cell set. Contract validation checks the source column and barcode overlap before the
analysis proceeds. See [advanced configuration](configuration.md#6-advanced-reproducibility-and-reusability).

## 6. Postprocessing and report

After annotating your cells/barcodes, you may run postprocessing rules:

```bash
snakemake run_postprocessing --sdm conda  # annotate → integrate → DE → neighbourhood → reports
```

This command will produce the integrated object, a per-sample report `results/integrated_samples/sample_report.pdf` to share with collaborators, pseudobulk DE (*optional*, but recommended), neighbourhood/co-occurrence and optional subcompartment re-clustering (for cell subtyping exploration).

## 7. Differential expression

SpaceBlocks performs differential expression on raw counts summed into pseudobulk observations. Comparison variables come from the `core_samples.tsv` sample sheet and may define a single contrast (for example, treated versus untreated), combined categorical groups (for example, phenotype and treatment), or a multilevel variable (for example, cancer stage).
Importantly, combined groups are compared directly, since SpaceBlocks does not currently fit statistical interaction terms. For multilevel variables, optional one-vs-rest tests identify level-specific biomarkers, while the omnibus LRT gives a holistic test across levels. Both require at least three eligible levels after exclusion, replication, and pairing filters; configured pairwise contrasts still run when either is skipped.

Independent samples should omit `paired_by`. When comparison levels are matched within the same patient or spatial sample, `paired_by` identifies that repeated unit. For each explicit pairwise contrast, only complete pairs are tested. See the [pseudobulk experimental-design configuration](configuration.md#pseudobulk-experimental-designs) for independent, combined-group, multilevel, and paired examples.

## 8. Explore genes and signatures

Once the integrated file is available, you may explore any gene or gene set of interest (generated in the results, or completely external).

List the genes / gene sets in the `gene_exploration.queries` TSV; integrated and per-sample plots share one expression scale for easy comparison.

```bash
snakemake explore_genes --sdm conda       # AUCell scoring + spatial/dotplot composites
```

## 9. Reproducibility

### 9.1 Recommendations at a glance

- **Colours** are consistent across output plots — set them once in `config.yaml`
  ([color scale configuration](configuration.md#4-color-scale-customization)).
- **Reproducibility**: to close the UMAP/Leiden non-determinism gap across machines, share
  pre-computed clusters/annotations or a pre-built contract h5ad (`mode: decoupled`) — see
  [reproducibility configuration](configuration.md#6-advanced-reproducibility-and-reusability).
- Make sure to back up the intermediate files (see section 9.2) and to [save a locked copy of your environment versions](environments.md).

### 9.2 Files to ensure reproducibility

To reproduce your full run, make sure to back up:
1. The design files (`config/core_samples.tsv`; `visiumhd_samples.csv`).
2. The GeoJSON files with the manual region annotations.
3. The metadata with the leiden resolutions, generated during preprocessing.
4. The annotated metadata, generated during the postprocessing. Optionally, if using manual annotation, the TSV mapping cluster-to-cell type equivalences.
5. The environment versions (recommended for publication).

## 10. Example runs

You may check our [public data end-to-end example runs](demos.md) to see how a full run looks.

Be mindful that wall time and memory usage depend strongly on dataset size and enabled analyses. Note that BANKSY spatial-niche construction is usually the most memory-intensive processing step.
