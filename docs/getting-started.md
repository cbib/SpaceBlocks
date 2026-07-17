# Get started: a full run

This page walks through a complete SpaceBlocks run and the recommendations to make it smooth. For the exhaustive key reference see [Configuration](configuration.md); for what each rule does, [Rule reference](rules.md).

If you want to see an example run of public data, refer to #TODO Visium HD example and Xenium 5K example.

> [!TIP]
> A run is defined almost entirely by two things: your **experimental design**
> (`core_samples.tsv`) and your **cell-type annotations**. Everything else has sensible
> defaults. Keep those two under version control and your run will be reproducible end-to-end.

##TODO INCLUDE NUMBERED INDEX

## 0. Prerequisites

- [Snakemake](https://snakemake.readthedocs.io) ≥ 8 and Conda/Mamba.
- For the **Visium HD** head only: an external [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger) ≥ 4.0.1.
- [QuPath](https://qupath.github.io/) (desktop), recommended to include manual region annotations for the downstream analyses (see [QuPath annotation tutorial](rules.md#qupath-annotation-tutorial)).
- Configure your Snakemake profile in `--profile profiles/default`. An example for HPCs with slurm scheduler is provided.

We recommend to always run snakemake with `--sdm conda`, so each rule gets its pinned environment.

## 1. Configure

In the `config/config.yaml` pick a `mode` and fill two files:

1. **`config/core_samples.tsv`** — one row per sample; extra columns (e.g. `patient`, `condition`) become design metadata, and can be used for downstream analyses and visualization.
2. **`config/config.yaml`** — set `mode`, the input and output paths, and (for Visium HD) the head files. See [Configuration](configuration.md) for full details.
3. If you are running `mode : visiumhd`, you will need to fill 

Before launching SpaceBlocks, we recommend testing your configuration using a dry-run:

```bash
snakemake -n --sdm conda      # dry-run: builds the DAG and validates the config
```
This will tell you if the config file has any incorrect input types (i.e. characters when expecting integers), and will build the environments for the different SpaceBlocks steps.


## 2. Region annotation (*optional*, but **recommended**)

When using any of the HeadBlocks, you may generate the images that will be included in the h5ad AnnData objects:

```bash
snakemake qupath_images --sdm conda
```
You may then manually annotate spatial regions in QuPath, export GeoJSONs and make them available for SpaceBlocks via `geojson_path` in `config/config.yaml`.

Then follow the [QuPath annotation tutorial](rules.md#qupath-annotation-tutorial), saving each file as `{sample}_tissue_hires_image.geojson` (Visium HD) or `{sample}_morphology.geojson` (Xenium 5K) into `config["geojson_path"]`. This step is optional —  Without it, regions are `Unlabeled`.

### 2.1. Region annotation in uncoupled mode

If you are using `mode : "uncoupled"`... #TODO

> [!WARNING]
> Make sure that the pixel size and coordinates of the images used to annotate your spatial regions are the same. This is automatically managed inside `qupath_images`, but may vary if you have a different output.

## 3. Preprocess, then look before you filter

SpaceBlocks provide a rule to screen how different QC filters (under `qc_sweep:` -> `min_genes`, `min_counts`, `max_counts`, `max_pct_mito`) would affect each of your samples.

You can combine this screening with an external annotated reference (`ingest_ref` and `ingest_ref_label_key` in the config) to observe how many, where in the tissue and which (automatically annotated) cell type would be removed when appliyin a threshold. To run it, use:

```bash
snakemake qc_sweep_all      --sdm conda   # OPTIONAL: see where each candidate threshold cuts
```

We recommend to run `qc_sweep` **first** to choose QC thresholds from your own data rather than the defaults.
You can next use the same thresholds for all samples or use a `sample_qc_thresholds` TSV if samples need different cut-offs.

After defining your general (*default*) or sample-specific thresholds, you should specify the preprocessing parameters (under `analysis:` in the config).

> [!TIP]
> We **strongly recommend** to include some literature predefined markers to screen for cell types present in your samples (defined via `snakemake_cell_markers`).

Next, you may run the preprocessing CoreBlock:

```bash
snakemake run_preprocessing --sdm conda   # validate → QC → normalise → cluster
```
The output (see [SpaceBlocks outputs](outputs.md) for details) will show spatial plots and UMAPs containing clustering metrics (clustree, Silhouette scores) for each resolution in the range, predefined markers (*optional*), QC plots to diagnose potentially noisy clusters. Then, you may proceed to the annotation step.

## 4. Annotate cell types

Annotation is the most limiting step in every high-throughput single cell pipeline. It conditions downstream analysis and, thus, it is a breakpoint between CoreBlocks within the SpaceBlocks workflow.

Through the development of SpaceBlocks, we have tested several automatic annotation tools in Visium HD data. The results in our own data had very limited accuracy and, thus, manual annotation is the default option. 

### 4.1 Manual annotation

SpaceBlocks automatically generates a cluster-to-cell type annotation TSV template within preprocessing. In this TSV, rows are cluster numbers and columns are sample names and Leiden resolution (i.e. "Sample1_res0.6").

The resolution is automatically read by the `annotation` rule, so in the TSV you need to update:
1. The equivalence between the cluster number and the cell type present.
2. The resolution desired to perform annotation in the column names.

After filling the annotation TSV, point `config["cluster_annotations"]` at it. This will unlock postprocessing rules.

### 4.2 External annotation

You may annotate your samples outside of SpaceBlocks, and continue the downstream analyses by providing annotated TSVs for each of them, see [Configuration §6](configuration.md#6-advanced-reproducibility-and-reusability).

## 5. Postprocessing and report

After annotating your cells/barcodes, you may run postprocessing rules:

```bash
snakemake run_postprocessing --sdm conda  # annotate → integrate → DE → neighbourhood → reports
```

This command will produce the integrated object, a per-sample PDF `sample_report` to share with collaborators, pseudobulk DE (*optional*, but recommended), neighbourhood/co-occurrence and optional subcompartment re-clustering (for cell subtyping exploration).

## 6. Explore genes and signatures

Once the integrated file is available, you may explore any gene or gene set of interest (generated in the results, or completely external).

List the genes / gene sets in the `gene_exploration.queries` TSV; integrated and per-sample plots share one expression scale for easy comparison.

```bash
snakemake explore_genes --sdm conda       # AUCell scoring + spatial/dotplot composites
```

## 7. Reproducibility

### 7.1 Recommendations at a glance

- **Colours** are consistent across every plot — set them once in `config.yaml`
  ([Configuration §4](configuration.md#4-color-scale-customization)).
- **Reproducibility**: to close the UMAP/Leiden non-determinism gap across machines, share
  pre-computed clusters/annotations or a pre-built contract (`mode: decoupled`) — see
  [Configuration §6](configuration.md#6-advanced-reproducibility-and-reusability).
- **Reproduce our published run** end-to-end on public data: [Reproduce the paper run](reproduction.md).

### 7.2 Files to ensure reproducibility

To reproduce your full run, make sure to back up:
1. The design files (`config/core_samples.tsv`; `samples_vhd.csv`).
2. The GeoJSON files with the manual region annotations.
3. The metadata with the leiden resolutions, generated during preprocessing.
4. The annotated metadata, generated during the postprocessing. Optionally, if using manual annotation, the TSV mapping cluster-to-cell type equivalences.
5. The environment versions (recommended for publication).
