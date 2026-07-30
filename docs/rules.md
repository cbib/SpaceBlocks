# SpaceBlocks — rule reference

Here you can find one-paragraph explanations for every rule, grouped by phase.

Inputs/outputs are summarised; the `.smk` files and `config["resources"]` are the exact source of truth.

!!! note "Optional rules"
    Rules marked *(optional)* only run when their config toggle / input is present.

---

## HeadBlocks

HeadBlocks are *optional* technology-specific modules, selected in `config["mode"]`. Each
one produces the standardized **contract h5ad** the CoreBlocks start from; the rules within
a headBlock carry a 3-letter technology suffix (`_vhd`, `_x5k`) so the organization stays
easy to follow.

The pipeline can also run in `mode: decoupled`, without any headBlock — the CoreBlocks then
consume pre-existing contract h5ads directly.

!!! tip "Annotate your regions first"
    Whichever headBlock you use, run the `qupath_images` target **first**: it produces the
    per-sample annotation images, which you annotate in QuPath and export as GeoJSON *before*
    launching the rest of the run. The contract builders (`prepare_input_*`) pick those
    GeoJSONs up automatically (and fall back to `Unlabeled` regions if none are present).

### Visium HD

Included only when `mode: visiumhd`. Produces the standardized contract h5ad.

!!! warning "Space Ranger installation requirement"
    This headBlock requires an external installation of [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/downloads) version ≥ 4.0.1 (proprietary 10x Genomics software), with its path set in `config["spaceranger"]`.

#### Visium HD: `spaceranger_count_vhd`
Runs Space Ranger (SR) `count` on each sample from the fastq files. An example input sheet can be found in `config/visiumhd_samples.csv`.
[Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/latest) is proprietary 10x Genomics software (path in `config["spaceranger"]`), so it must be installed externally.
Output is a `.done` sentinel; the real SR files are tracked by the rules that consume them, giving Snakemake correct provenance without Martian conflicts.

#### Visium HD: `generate_qupath_vhd`
A cheap file copy: finds `tissue_hires_image.png` in the version-varying Space Ranger tree and copies it to `Samples/{sample}/QuPath_image/`.

This image is meant to be annotated in QuPath to generate a GeoJSON with the region-annotation polygons — see the [QuPath annotation tutorial](qupath-tutorial.md).
The GeoJSON files follow the naming convention `{sample}_tissue_hires_image.geojson` and are read from the folder in `config["geojson_path"]`.

#### Visium HD: `prepare_input_vhd`
Takes the `spaceranger_count_vhd` output and the QuPath GeoJSON region annotations, and prepares the contract h5ad for the CoreBlocks.
Discovers the SR files, reads the raw count matrix, computes per-cell centroids (mean position of the 2 µm bins making up each segmented cell), embeds the hires image + scalefactors into `uns["spatial"]`,
joins the QuPath GeoJSON into `obs["region_annotation"]`, and writes the **unfiltered contract h5ad** (raw counts, no filtering/normalisation). This is the file every core rule starts from.

### Xenium 5K

Included only when `mode: xenium5k`. Produces the *same* standardized contract h5ad as the Visium HD head, so the CoreBlocks run identically. The head only needs the SpatialData stack (`envs/xenium5k.yaml`).

#### Xenium 5K: `convert_zarr_x5k`
The heavy I/O step (once per sample): converts a Xenium output bundle to a [SpatialData](https://spatialdata.scverse.org/) zarr store via `spatialdata_io`, which becomes the canonical object holding the morphology image, cell masks, boundaries, and the raw count table.
Analogous to `spaceranger_count_vhd`. The tracked output is a `{sample}.zarr.done` marker that the script writes **only after** the zarr is fully written, so an interrupted conversion leaves no marker and is correctly re-run; the downstream rules derive the zarr path by stripping `.done`.

#### Xenium 5K: `generate_qupath_x5k`
Reads the `morphology_focus` channels from the zarr at a configurable pyramid level and composites them into an RGB TIFF for QuPath, plus a scale-factor JSON for QuPath-pixel → micron coordinate mapping.
Analogous to `generate_qupath_vhd`. Channel colours are assigned by **name** (nuclear/DAPI → blue, membrane/boundary → red, RNA → green, cytoplasm → amber), with a positional-palette fallback for unrecognised channels and any channel count, and the mapping is logged each run so it can be verified. The GeoJSON files follow the naming convention `{sample}_morphology.geojson` (read from `config["geojson_path"]`).

#### Xenium 5K: `prepare_input_x5k`
Takes the zarr and the QuPath GeoJSON, and writes the **unfiltered contract h5ad** for the CoreBlocks.
Extracts the raw count table, embeds a channel-agnostic greyscale morphology image + scalefactors into `uns["spatial"]`, joins the GeoJSON into `obs["region_annotation"]`, and writes raw counts with **no filtering/normalisation/clustering** — that is deliberately the CoreBlock's job, so a Xenium run and a Visium HD run reach the core identically. Analogous to `prepare_input_vhd`.

---

## CoreBlocks

The CoreBlocks are the central (core) modules of SpaceBlocks.

They perform the `preprocessing`, `postprocessing` and `exploration` steps for single-cell resolution spatial transcriptomics data, starting from a **standardized ("contract") `AnnData` h5ad file**.

### CoreBlock 1 — preprocessing

The first CoreBlock includes rules from data preparation to clustering.

#### `validate_input`
The **most important rule** of the pipeline: performs contract-format validation.

Confirms that the contract for each sample in `config/core_samples.tsv` is well-formed (loads as AnnData; >0 cells/genes; non-negative integer raw counts; finite `obsm["spatial"]`; `sample` in obs; unique names; well-formed `uns["spatial"]` if present).
It supports two modes: requiring region annotations (hard) or not (soft). Hard failures prevent the DAG from advancing. Soft issues (missing region annotation, no mito genes, no image) are logged into `validation/input_validation.json` and are not fatal.

#### `qc_sweep` *(optional; target `qc_sweep_all`)*
Optional rule that provides diagnostic information about QC before committing to downstream analyses. It never filters, clusters, or writes an h5ad.
The output shows *where on the tissue* each candidate QC threshold would remove cells, before any filter is committed.
It includes joint scatter, violin, and spatial **plots** plus a **summary TSV** showing how the cut-offs would affect each sample.
Automatic annotation from an externally annotated dataset (see `ingest_ref`) is optional, shown side-by-side with the spatial output.

#### `preprocess_umap`
The **core preprocessing step**, technology-agnostic.

Reads the validated contract, applies QC filters (shared or per sample) — or, with external annotation, masks to the externally-annotated cells; normalises + log1p (keeping a raw-counts layer), runs PCA → neighbours → UMAP → multi-resolution Leiden, and stamps core-sheet metadata into `obs`.
Writes `adata_{sample}.h5ad`, `metadata_{sample}.tsv`, and a per-sample `{sample}_report.tsv` (cells before/after, per-feature min/max, applied cut-offs).

#### `leiden_analysis`
Clustering diagnosis and marker visualization before (manual) annotation.

Parallelizes Leiden cluster plots per sample and resolution (computed in `preprocess_umap`). Produces individual and merged cluster plots (UMAP + spatial) for each Leiden resolution, silhouette scores, clustree plots, and violin plots for the common QC metrics.

#### `generate_annotation_template`
Writes an empty TSV template to facilitate writing the `annotate_cells` input, supporting in-pipeline manual annotation.
You fill in cell-type labels and point `config["cluster_annotations"]` at it to unlock annotation and the annotation-dependent downstream rules.

#### `ingest_ref` *(optional)*
Automatic annotation via `scanpy.ingest()` from an annotated scRNA-seq reference. Facilitates manual annotation by comparing each sample against an external reference.

#### `spatial_niches` *(optional)*
Runs BANKSY across concatenated samples to identify shared spatial niches.
Produces per-sample niche TSVs (barcode → niche), a concatenated h5ad, and spatial plots for each niche.

---

### CoreBlock 2 — postprocessing

This second CoreBlock includes rules from sample annotation to integration and differential expression.

#### `annotate_cells`
Adds cell-type annotations from a manually curated TSV (cluster → cell type), or from an external annotation TSV (barcode → cell type).
Writes `adata_{sample}_annotated.h5ad` plus annotation plots and composition barplots.

#### `integrate_samples`
Concatenates all annotated samples, runs Harmony batch correction on the specified variable (`integration.integrate_key`), geosketch subsampling, and composition barplots.

Writes `concatenated.h5ad`, `harmony_integrated.h5ad`, and `sketched.h5ad`. This is the input to pseudobulk, subclustering (CoreBlock 2), and gene exploration (CoreBlock 3).

#### `pseudobulk_aggregate`
Builds pseudobulk count matrices via [decoupler](https://decoupler.readthedocs.io/en/latest/) per `region`, `cell type` and (optionally) `niche` from the integrated object.

Produces diagnostic QC and PCA plots coloured by sample and by each extra-annotation column.

#### `pseudobulk_de` *(optional; `analysis.run_pseudobulk_de: true`)*
Identifies deregulated genes between conditions and annotated regions via differential expression with DESeq2 (R package) on the pseudobulk matrices.

Differential expression runs in three ways:
- pairwise Wald contrasts between level pairs (region A vs region B),
- each level against the rest,
- and, if there are more than 2 levels, a likelihood-ratio test (LRT) with DEGpatterns for gene-group clustering.

Produces TSV tables, volcano plots, and metadata-annotated heatmaps between conditions.

#### `neighbourhood_analysis`
Runs squidpy neighbourhood enrichment per sample and per `region`, plus per-region cell-type co-occurrence (a full co-occurrence figure per region and compact pairwise heatmaps).

#### `subcluster` *(optional)*
Re-clustering of the specified subcompartment(s) (a config-specified subset of cell types) for cell-subtype exploration.

Generates integrated (Harmony) and unintegrated results, with UMAPs coloured by `region`, spatial niche, and external metadata columns, marker tables, and composition barplots.

#### `sample_report`
Assembles a compact, multi-page PDF with graphical outputs per sample, ideal to share with internal and external collaborators.

Each page contains the most important data-descriptive information: cluster / annotation / region / niche UMAPs and matching spatial maps, bar plots for absolute and relative composition, and a dot plot with the top-10 expressed markers per cell type.

---

### CoreBlock 3 — exploration

This third CoreBlock allows the informative exploration of curated genes and gene sets (obtained in the DE step or from external sources) in individual and merged samples.

Integrated exploration scales are standardized to ease interpretation of the results.

#### `explore_genes_integrated`
Exploration of the integrated samples.

Loads the Harmony-integrated object once and computes AUCell scores for the input gene or gene set. It then writes shared expression ranges and produces integrated UMAP/dotplot PNGs per gene or gene set.

#### `explore_genes_sample`
Exploration of individual samples.

Using the shared expression ranges from the integrated exploration, produces per-sample spatial composites (expression + cell type + region) and stacked dotplot composites for each query entry.

---

## Resources & retries

No rule uses `localrule`, so nothing runs on the scheduler's head node.

Every compute rule draws `mem_mb` / `runtime` / `threads` from `config["resources"]` (with a `default` fallback).

Memory is wrapped so it **grows with the retry attempt** (`mem_mb = base × attempt`); combined with `retries: 3` in the SLURM profile, an OOM-killed job is automatically resubmitted with more RAM.
