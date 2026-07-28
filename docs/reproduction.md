# Reproduce the paper run — `mode: decoupled` demos

This page shows how to run SpaceBlocks end-to-end on **public data**, and documents how to prepare your own inputs for `mode: decoupled`.

> [!WARNING]
> The tiny case in `.tests/` exists only to exercise the DAG in CI, and to render the workflow's tube map on the Snakemake catalog.

- [Technical notes](#technical-notes)
- [Visium HD example (mouse brain)](#visium-hd-example-mouse-brain)
- [Xenium 5K example (human melanoma)](#xenium-5k-example-human-melanoma)
  - [Configure and run the examples](#configure-and-run-the-examples)
  - [What preconfigured files exercise](#what-preconfigured-files-exercise)
- [Preparing other inputs for decoupled mode](#preparing-inputs-for-decoupled-mode)


## Technical notes

To keep the tutorials simple and lightweight, we use 1 sample from a public dataset and divide it into **3 artificial samples, keeping only 500 HVGs**.

We provide preconfigured `config` files, pre-annotated regions (GeoJSON files), clusters and cluster-to-cell type equivalences for the example public datasets to ensure reproducibility.

The exact environment versions used during the generation of these tutorials can be found under `reproduction/lock.envs`.

Each of the examples below can run under 8Gb of RAM and 4 cores in <40 min under our `slurm` HPC system:

```
Our HPC specifications

7 standard compute nodes (HP Apollo 2000 Gen10+)
2 × AMD EPYC 7513 CPUs (32 cores, 2.8 GHz each)
256 GB RAM per node

1 high-memory compute node (HP Apollo 2000 Gen10+)
2 × AMD EPYC 7513 CPUs (32 cores, 2.8 GHz each)
1 TB RAM
```

## Visium HD example (mouse brain)

SpaceBlocks is built to analyse single-cell resolution Spatial Transcriptomics data, so Visium HD data needs to be preprocessed via [bin2cell](https://github.com/Teichlab/bin2cell), [ENACT](https://github.com/Sanofi-Public/enact-pipeline) or, as in the Visium HD HeadBlock, Space Ranger >= v4.0.1 (internally implementing StarDist segmentation).

> [!IMPORTANT]
> We use a **Space Ranger ≥ 4.0.1** Visium HD dataset (one that ships `segmented_outputs/`). The older Space Ranger 3.x Mouse Brain release has **no** segmentation, so `format_visiumhd.py` would > find no cell table.

For the example here presented, you may download the Visium HD dataset from the [10x Genomics web](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4), or via terminal using curl or wget. This dataset has been processed with Space Ranger v4.0.1:

```bash
# at the workflow dir, build the environment
cd SpaceBlocks/
conda env create --file=workflow/envs/visiumhd.yaml
# activate it
conda activate visiumhd
# move to the reproducibility dir and place the downloaded outs under data/<SAMPLE>/
# move to the reproducibility dir and download the three files under data/<SAMPLE>/
cd reproduction/visiumhd
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
```

We provide `demo_vhd.geojson` for this dataset, as an example QuPath export on the hires image (see [QuPath annotation tutorial](qupath-tutorial.md)).

## Xenium 5K example (human melanoma)

You may download the Xenium dataset from the [10x Genomics web](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin), or via terminal using curl or wget:

```bash
# at the workflow dir, build the environment
cd SpaceBlocks/
conda env create --file=workflow/envs/xenium5k.yaml
# activate it
conda activate xenium5k
# move to the reproduciblity dir
cd reproduction/xenium5k
# download the Xenium bundle into data/<SAMPLE>/
mkdir -p data && cd data
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Human_Skin_FFPE/Xenium_Prime_Human_Skin_FFPE_outs.zip
unzip Xenium_Prime_Human_Skin_FFPE_outs.zip -d Xenium_Prime_Human_Skin_FFPE
cd ..
# build the contracts
# use an interactive job with 16Gb and 4 cores if using a slurm cluster (salloc --ntasks=1 --cpus-per-task=4 --mem 16G -t 01:00:00)
python format_xenium.py              # -> contracts/<roi>.h5ad + core_samples.tsv
```

Notice that the Xenium GeoJSON is annotated in pixels, while cells are in microns. The `format_xenium.py` script detects the scale mismatch and rescales the polygons automatically.

Importantly, **`format_xenium.py` embeds a greyscale composite of the `morphology_focus` channels in `uns["spatial"]`**, so spatial plots are drawn over the tissue instead of on a bare scatter. Set `HIRES_LEVEL` in the script to trade resolution for file size.

## Configure and run the examples

A ready-to-run config ships next to each example: `reproduction/xenium5k/xenium5k_config.yaml` and `reproduction/visiumhd/visiumhd_config.yaml`. Optionally, you may use an external reference for the `ingest` rule, we provide external links to set in the example config files.

To use one, comment the default `configfile:` line in `workflow/Snakefile` and uncomment the matching example line (all are already present):

```python
# configfile: "config/config.yaml"
# configfile: reproduction/visiumhd/visiumhd_config.yaml      # alternative for visiumhd
configfile: "reproduction/xenium5k/xenium5k_config.yaml"
```

> [!WARNING]
> **Make sure to only have one uncommented config file in the Snakefile.**
> Snakemake allows the Snakefile to reference multiple config files, so having multiple config files would lead to unexpected results.

Then launch the core:

```bash
snakemake run_preprocessing qc_sweep_all --sdm conda
snakemake run_postprocessing             --sdm conda   # integration + DE across the 3 synthetic samples
```

### What preconfigured files exercise

#### Config files

The shipped `reproduction/xenium5k/xenium5k_config.yaml` is set up to demonstrate the optional rules too:

- **Subcompartment re-clustering.** `subcompartments` groups cell type labels and re-clusters them (`subcluster` rule, included in `run_posprocessing`).
- **Gene / signature exploration.** `gene_queries_demo.tsv` holds a small gene set and a gene to showcase the `run_exploration` CoreBlock capabilities.
- **Automatic annotation (optional).** `ingest_ref` and `ingest_ref_label_key` are empty by default, they allow the use of external references for automatic annotation (which may guide manual annotation). We provide examples in the config files to test this feature during `qc_sweep_all` (showing which assigned cell types are removed during QC) and `run_preprocessing`. Automatic annotation will also be shown in the `sample_report.pdf` during `run_postprocessing`.

#### TSV cluster annotation files

Manual assignation of Leiden clusters into cell types (via `config/cluster_annotations.tsv`) is the default requirement to run the postprocessing and exploration CoreBlocks in SpaceBlocks.

Example cell type annotation files are provided as `cluster_annotations_{technology}_demo.tsv` to show how the TSV should be filled.

## Preparing inputs for decoupled mode

The config `mode: decoupled` skips the HeadBlocks entirely and runs the CoreBlocks on contract h5ads you provide. The two scripts above are worked examples of how to produce them (with artificial samples).

To generate valid *contract* AnnData objects (saved as h5ads), follow the [contract specification](design.md#the-validation-object-contract); the rules in short are:

**1 — Shape.** Each h5ad must satisfy the [contract](design.md#the-validation-object-contract): `X` = raw integer counts, `obsm["spatial"]` = per-cell `(x, y)`, and `obs` carrying `sample` (and, ideally, `cell_id` and `region_annotation`). `cell_id` is **optional but strongly recommended, and must equal `obs_names`** to prevent unexpected behaviour — downstream joins key on it, and `validate_input` warns when it differs. An embedded `uns["spatial"]` image is optional, but also strongly recommended. `validate_input` checks all of this before anything runs.

**2 — Naming & layout.** One file per sample, laid out as:

```
<contract_dir>/<sample>.h5ad
```

where each `<sample>` is a row in `config/core_samples.tsv`. The sample names in the sheet, the folder names, and the `obs["sample"]` values must all match. The example scripts (see above) guarantee this. You can find the TSV examples in `reproduction/<technology>/contracts/core_samples.tsv`.

**3 — Config.** Point the run at them:

```yaml
mode: "decoupled"
contract_dir: "reproduction/visiumhd/contracts"   # holds <sample>.h5ad
core_samples: "reproduction/visiumhd/contracts/core_samples.tsv"
# plus post_processing_outdir, logdir, geojson_path (see Configuration)
```

Region annotation in decoupled mode is baked into the contract (as `region_annotation`) when you
build it — there is no `qupath_images` target and no GeoJSON join later.





