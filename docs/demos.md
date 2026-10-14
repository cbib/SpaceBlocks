# Demos — public-data runs per technology

This page shows how to run SpaceBlocks end-to-end on **public data** for each supported technology, and documents how to prepare your own inputs for `mode: decoupled`.

All demo config paths are relative to the repository root: place downloaded inputs under `demos/<platform>/data/` and run Snakemake from the repository root.

!!! warning "Repository's synthetic data"
    The tiny case in `.test/` exists only to exercise the DAG in CI, and to render the workflow's tube map on the Snakemake catalog.

## Technical notes

To keep the Visium HD and Xenium 5K tutorials simple and lightweight, we use 1 sample from a public dataset and divide it into **3 artificial samples, keeping only 500 highly-variable genes (HVGs)**. The Atera and MERSCOPE examples consume publicly available raw data.

We provide default `config` files, pre-annotated regions (GeoJSON files), clusters and cluster-to-cell type equivalences for the example public datasets to ensure reproducibility.

The exact environment versions used during the generation of these tutorials can be found under `demos/lock.envs`.

The small decoupled Visium HD and Xenium 5K examples run under 8 GB of RAM and four cores
in under 40 minutes on our Slurm system. The full HeadBlock examples have platform- and
dataset-dependent requirements:

```
Our HPC specifications

7 standard compute nodes (HP Apollo 2000 Gen10+)
2 × AMD EPYC 7513 CPUs (32 cores, 2.8 GHz each)
256 GB RAM per node

1 high-memory compute node (HP Apollo 2000 Gen10+)
2 × AMD EPYC 7513 CPUs (32 cores, 2.8 GHz each)
1 TB RAM
```

## Visium HD decoupled mode example (mouse brain)

SpaceBlocks is built to analyse single-cell resolution Spatial Transcriptomics data, so Visium HD data needs to be preprocessed via [bin2cell](https://github.com/Teichlab/bin2cell), [ENACT](https://github.com/Sanofi-Public/enact-pipeline) or, as in the Visium HD HeadBlock, Space Ranger >= v4.0.1 (internally implementing StarDist segmentation).

!!! important "Space Ranger version matters"
    We use a **Space Ranger ≥ 4.0.1** Visium HD dataset (one that ships `segmented_outputs/`). The older Space Ranger 3.x Mouse Brain release has **no StarDist segmentation**, so `format_visiumhd.py` would find no cell table.

For the example here presented, you may download the Visium HD dataset from the [10x Genomics web](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4), or via terminal using `curl` or `wget`. This dataset has been processed with Space Ranger v4.0.1.

!!! note "Visium HD demo dataset size and SpaceBlocks mode"
    This demo runs in `mode: decoupled` and the size of the dataset is ~15GB.

```bash
# from the repository root, build and activate the environment
conda env create --file=workflow/envs/visiumhd.yaml
conda activate visiumhd
# move to the demo directory and download the three files under data/<SAMPLE>/
cd demos/visiumhd
mkdir -p data/Visium_HD_Mouse_Brain && cd data/Visium_HD_Mouse_Brain
BASE=https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Mouse_Brain
curl -O $BASE/Visium_HD_Mouse_Brain_segmented_outputs.tar.gz
curl -O $BASE/Visium_HD_Mouse_Brain_barcode_mappings.parquet
curl -O $BASE/Visium_HD_Mouse_Brain_binned_outputs.tar.gz
tar -xzf Visium_HD_Mouse_Brain_segmented_outputs.tar.gz    # -> segmented_outputs/
tar -xzf Visium_HD_Mouse_Brain_binned_outputs.tar.gz       # -> binned_outputs/
cd ../..
# build the contracts for the CoreBlocks
# use an interactive job with 16Gb and 4 cores on a slurm cluster (salloc --ntasks=1 --cpus-per-task=4 --mem 16G -t 01:00:00)
python format_visiumhd.py            # -> contracts/<sample>.h5ad + core_samples.tsv
cd ../..                             # return to the repository root
```

We provide `demo_vhd.geojson` for this dataset, as an example QuPath export on the hires image (see [QuPath annotation tutorial](qupath-tutorial.md)).

## Xenium 5K decoupled mode example (human melanoma)

You may download the Xenium dataset from the [10x Genomics web](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin), or via terminal using curl or wget.

!!! note "Xenium 5K demo dataset size and SpaceBlocks mode"
    This demo runs in `mode: decoupled` and the size of the dataset is ~11GB.

```bash
# from the repository root, build and activate the environment
conda env create --file=workflow/envs/xenium5k.yaml
conda activate xenium5k
# move to the demo directory
cd demos/xenium5k
# download the Xenium bundle into data/<SAMPLE>/
mkdir -p data && cd data
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Human_Skin_FFPE/Xenium_Prime_Human_Skin_FFPE_outs.zip
unzip Xenium_Prime_Human_Skin_FFPE_outs.zip -d Xenium_Prime_Human_Skin_FFPE
cd ..
# build the contracts
# use an interactive job with 16Gb and 4 cores if using a slurm cluster (salloc --ntasks=1 --cpus-per-task=4 --mem 16G -t 01:00:00)
python format_xenium.py              # -> contracts/<roi>.h5ad + core_samples.tsv
cd ../..                             # return to the repository root
```

Notice that the Xenium GeoJSON is annotated in pixels, while cells are in microns. The `format_xenium.py` script detects the scale mismatch and rescales the polygons automatically.

Importantly, **`format_xenium.py` embeds a greyscale composite of the `morphology_focus` channels in `uns["spatial"]`**, so spatial plots are drawn over the tissue instead of on a bare scatter. Set `HIRES_LEVEL` in the script to trade resolution for file size.

## Atera HeadBlock example (human breast cancer, alpha)

!!! warning "Alpha support"
    Atera is now shipping, but the dataset below remains a public **preview** generated with
    development software and a pre-release output format. This is a HeadBlock compatibility
    test rather than a polished production-data demo.

Unlike the two examples above, this one runs a HeadBlock rather than `mode: decoupled`: the contract h5ad is built by `prepare_input_ate` from the raw bundle, so it exercises the whole head chain. The [dataset](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer) is 18,028 genes over 170,057 cells (CC BY 4.0).

!!! note "Atera demo dataset size and SpaceBlocks mode"
    This demo runs in `mode: atera` and the size of the dataset is ~190GB, and the total space used after running SpaceBlocks is close to 300GB.

```bash
# from the repository root, build and activate the environment
conda env create --file=workflow/envs/atera.yaml
conda activate atera
# download into the configured demo data directory (allow ~300 GB free)
mkdir -p demos/atera/data && cd demos/atera/data
BASE=https://s3-us-west-2.amazonaws.com/10x.files/samples/atera/dev/WTA_Preview_FFPE_Breast_Cancer
curl -O $BASE/WTA_Preview_FFPE_Breast_Cancer_outs.zip
unzip WTA_Preview_FFPE_Breast_Cancer_outs.zip -d outs
# supplemental annotation files used by the demo
SUP=https://cf.10xgenomics.com/samples/atera/dev/WTA_Preview_FFPE_Breast_Cancer
curl -O $SUP/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv
curl -O $SUP/WTA_Preview_FFPE_Breast_Cancer_gene_groups.csv
# registered H&E and alignment (clear the three he_* config keys if omitted)
curl -O $SUP/WTA_Preview_FFPE_Breast_Cancer_he_image.ome.tif
curl -O $SUP/WTA_Preview_FFPE_Breast_Cancer_he_alignment.csv
curl -O $SUP/WTA_Preview_FFPE_Breast_Cancer_keypoints.csv
cd ../../..                         # return to the repository root
```

The shipped config already points to these paths. It annotates cells directly from 10x's per-cell calls (`annotation_types: [external_annotation]`), so the run needs no manual cell-type annotation step. Convert the supplemental cell-group file with the helper that ships next to the config:

```bash
# 10x's cell calls -> external-annotation metadata
python demos/atera/format_atera.py metadata \
    demos/atera/data/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv \
    --sample WTA_Preview_FFPE_Breast_Cancer --outdir demos/atera/external_metadata

snakemake qupath_images qc_sweep_all --sdm conda   # optional regions: annotate a QuPath TIFF and export the GeoJSON
snakemake run_preprocessing          --sdm conda
snakemake run_postprocessing         --sdm conda
```

The helper also prints the `annotation_colors` block using 10x's own display colours, already pasted into the shipped config.

Region annotation is optional here (the demo annotates cells from 10x's own calls). If you do want regions, `qupath_images` exports two annotatable TIFFs per sample, and `prepare_input_ate` reads whichever GeoJSON you export **named after the image stem** (`{sample}_<image>.geojson`), preferring them in this order:

- `{sample}_he_background.tiff` — the registered H&E resampled onto the morphology grid (also the embedded contract background). **Recommended**: it sits on the cells' coordinate grid, so regions map back with a single scalar (no affine). Export as `{sample}_he_background.geojson`.
- `{sample}_morphology.tiff` — the fluorescence composite (the fallback when no H&E is configured). Export as `{sample}_morphology.geojson`.

Importantly, the registered H&E is not the raw whole-slide scan, which covers a second tissue section.

The bundle holds a single sample, so `run_postprocessing` is not meaningful/complete here because integration, pseudobulk DE and composition comparisons all need replicates. The shipped config therefore trims those options; see `demos/xenium5k/xenium5k_config.yaml` for a config that exercises them.

## MERSCOPE HeadBlock example (mouse brain)

Like the Atera demo, the MERSCOPE demo runs its HeadBlock. The contract h5ad is built by `prepare_input_mer` straight from a raw Vizgen region, so it exercises the whole MERSCOPE head chain (`generate_qupath_mer` → `prepare_input_mer`). The head reads the two Vizgen cell CSVs directly and reuses the `xenium5k.yaml` environment (no intermediate Zarr and no new environment).

!!! note "MERSCOPE demo dataset size and SpaceBlocks mode"
    This demo runs in `mode: merscope` and the size of the dataset is ~25GB.


We use the open-access **Vizgen MERFISH Mouse Brain Receptor Map** (Vizgen Data Release V1.0), Slice 1 / Replicate 1 — 483 genes, ~78k cells, a full coronal slice at single-cell resolution. Vizgen distributes this dataset through its **Data Release Program** landing page rather than an open CDN, so there is no single `curl` URL: request/download the files from [info.vizgen.com/mouse-brain-map](https://info.vizgen.com/mouse-brain-map) (the underlying Google Cloud bucket is requester-pays)

```bash
# from SpaceBlocks/, set the data Google Cloud dirs
BASE="gs://public-datasets-vizgen-merfish/datasets/mouse_brain_map/BrainReceptorShowcase/Slice1/Replicate1"
DEST=demos/merscope/data/BrainReceptorShowcase_S1R1
mkdir -p "$DEST/images"
# download the data for the demo using gcloud
gcloud storage ls "$BASE/images/"
gcloud storage cp "$BASE/cell_by_gene_S1R1.csv"   "$DEST/"
gcloud storage cp "$BASE/cell_metadata_S1R1.csv"  "$DEST/"
gcloud storage cp "$BASE/images/micron_to_mosaic_pixel_transform.csv" "$DEST/images/"
gcloud storage cp "$BASE/images/mosaic_DAPI_z3.tif"  "$DEST/images/"
gcloud storage cp "$BASE/images/mosaic_PolyT_z3.tif" "$DEST/images/"
```

After the download, the data should look like this:

```text
demos/merscope/data/BrainReceptorShowcase_S1R1/
├── datasets_..._cell_by_gene_S1R1.csv        # raw transcripts-per-cell matrix
├── datasets_..._cell_metadata_S1R1.csv       # centroids (center_x/y, µm) + QC
└── images/
    ├── datasets_..._micron_to_mosaic_pixel_transform.csv
    ├── ..._mosaic_DAPI_z3.tif                 # nuclear stain  (always present)
    └── ..._mosaic_PolyT_z3.tif               # total-RNA stain (always present)
```

The mouse brain release predates the Cell Boundary Stain Kit, so it ships **DAPI + PolyT only**. The MERSCOPE HeadBlock composites whatever channels it finds and skips absent ones, so no config change is needed.

If your download names its mosaics with a z-plane other than `z3`, set `merscope.z_index` to a plane that exists (the head errors clearly if the requested channel/z is missing).

```bash
# from the repository root, build and activate the environment (shared with Xenium)
conda env create --file=workflow/envs/xenium5k.yaml
conda activate xenium5k

# the shipped config already points to the region above; select it in workflow/Snakefile
# (see "Configure and run" below), then run:
snakemake qc_sweep_all qupath_images       --sdm conda   # OPTIONAL qc_sweep and regions: annotate the composite TIFF,
                                                         # export demos/merscope/geojson/<sample>_morphology.geojson
snakemake run_preprocessing                --sdm conda   # builds the contract via the head, then the CoreBlock
```

`generate_qupath_mer` writes `<sample>_morphology.tiff` (the annotation image), a grey background, and a `<sample>_morphology_scalefactors.json` that records the mosaic-affine mapping **once**; `prepare_input_mer` reads it back so cell centroids, the embedded background and any QuPath regions all share one grid (see [design](design.md#the-validation-object-contract) for the contract shape).

Region annotation is optional. Without a GeoJSON, `region_annotation` is `Unlabeled`.

The shipped config holds a **single** region, so postprocessing (integration, pseudobulk DE, composition comparisons) is not meaningful and has been trimmed from the `merscope_config.yaml`, exactly as in the Atera demo. To exercise the full core, add Slice 1 Replicate 2 and Replicate 3 as extra rows in `demos/merscope/core_samples.tsv` (they are genuine biological replicates) and re-enable the postprocessing options (see `demos/xenium5k/xenium5k_config.yaml` for a config that uses them).

## Configure and run the examples

A config ships next to each example: `demos/xenium5k/xenium5k_config.yaml`, `demos/visiumhd/visiumhd_config.yaml`, `demos/atera/atera_config.yaml` and `demos/merscope/merscope_config.yaml`. Their default paths match the layouts shown above when Snakemake is run from the repository root.

The Xenium 5K config also demonstrates optional automatic annotation. It expects the annotated reference named by `ingest_ref` under `demos/xenium5k/data/`; either provide that file or disable the feature by setting `qc_sweep.ingest_enabled: false` and `ingest_ref: ''`.

To use one, comment the default `configfile:` line in `workflow/Snakefile` and uncomment the matching example line (all are already present):

```python
# configfile: "config/config.yaml"
# configfile: "demos/visiumhd/visiumhd_config.yaml"    # alternative for visiumhd
# configfile: "demos/atera/atera_config.yaml"          # Atera preview-data HeadBlock test
# configfile: "demos/merscope/merscope_config.yaml"    # MERSCOPE HeadBlock test
configfile: "demos/xenium5k/xenium5k_config.yaml"
```

!!! warning "DO NOT INCLUDE several config files"
    **Make sure only one `configfile:` line is uncommented in the Snakefile.** Snakemake can merge multiple config files, which may produce unexpected results here.

Then launch the core:

```bash
snakemake run_preprocessing qc_sweep_all --sdm conda
snakemake run_postprocessing             --sdm conda   # integration + DE across the 3 synthetic samples
```

### What preconfigured files exercise

#### Config files

The shipped `demos/xenium5k/xenium5k_config.yaml` is set up to demonstrate the optional rules too:

- **Subcompartment re-clustering.** `subcompartments` groups cell type labels and re-clusters them (`subcluster` rule, included in `run_postprocessing`).
- **Gene / signature exploration.** `gene_queries_demo.tsv` holds a small gene set and a gene to showcase the `run_exploration` CoreBlock capabilities.
- **Automatic annotation (optional).** The Xenium 5K config enables this example and points `ingest_ref` to `demos/xenium5k/data/bcc_and_normal_Garnier_PNAS_2024.h5ad`; `ingest_ref_label_key` identifies the reference's cell-type column. The assigned cell types appear during `qc_sweep_all`, `run_preprocessing`, and in `sample_report.pdf` during `run_postprocessing`. To skip it, set `qc_sweep.ingest_enabled: false` and `ingest_ref: ''`.

#### TSV cluster annotation files

Manual assignation of Leiden clusters into cell types (via `config/cluster_annotations.tsv`) is the default requirement to run the postprocessing and exploration CoreBlocks in SpaceBlocks.

Example cell type annotation files are provided as `cluster_annotations_{technology}_demo.tsv` to show how the TSV should be filled.

## Preparing inputs for decoupled mode

The config `mode: decoupled` skips the HeadBlocks entirely and runs the CoreBlocks on contract h5ads you provide. The two scripts above are worked examples of how to produce them (with artificial samples).

To generate valid *contract* AnnData objects (saved as h5ads), follow the [contract specification](design.md#the-validation-object-contract); the rules in short are:

**1 — Shape.** Each h5ad must satisfy the [contract](design.md#the-validation-object-contract): `X` = raw integer counts, `obsm["spatial"]` = per-cell `(x, y)`, and `obs` carrying `sample` (and, ideally, `cell_id` and `region_annotation`). `cell_id` is **optional but strongly recommended**; when present, it must equal `obs_names` to prevent unexpected behaviour. Downstream joins key on it, and `validate_input` warns when it differs. An embedded `uns["spatial"]` image is optional, but also strongly recommended. `validate_input` checks all of this before anything runs.

**2 — Naming & layout.** One file per sample, laid out as:

```
<contract_dir>/<sample>.h5ad
```

where each `<sample>` is a row in `config/core_samples.tsv`. The sample names in the sheet, the `.h5ad` file names, and the `obs["sample"]` values must all match. The example scripts above guarantee this. You can find the generated TSV examples in `demos/<technology>/contracts/core_samples.tsv`.

**3 — Config.** Point the run at them:

```yaml
mode: "decoupled"
contract_dir: "demos/visiumhd/contracts"   # holds <sample>.h5ad
core_samples: "demos/visiumhd/contracts/core_samples.tsv"
# plus post_processing_outdir and logdir (see Configuration)
```

Region annotation in decoupled mode is baked into the contract (as `region_annotation`) when you build it.

External cell-type labels may likewise be stored as an `obs` column in each contract. Set
`external_annotation.enabled: true` and point `external_annotation.column` at that column; a
separate `precomputed_metadata_dir` is optional for this decoupled-only input path. Contract
`obs_names` remain the join key, and validation checks that the configured labels are present.
