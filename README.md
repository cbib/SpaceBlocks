# SpaceBlocks: streamlined, reproducible, spatial transcriptomics analyses

SpaceBlocks is a technology-agnostic Snakemake workflow to **semi-automatically analyse single-cell-resolution spatial transcriptomics (ST).**

SpaceBlocks allows **easy, customisable, shareable and fully reproducible** quality control, clustering, annotation, integration, differential expression, and spatial neighbourhood/niche analysis.

The architechture is meant to allow the analysis of public and *in house* generated ST datasets, with a modular and long-term maintainable design.


## Overview

The workflow is divided into (optional) technology-specific **HeadBlocks** and common **CoreBlocks** that streamline pre-, post-processing, and informative exploration of results.

<p align="center"><img src="docs/img/overview.svg" alt="SpaceBlocks workflow overview" width="760"></p>
 
A full SpaceBlocks run takes three inputs:

1. **ST formatted AnnData objects**. Generated from the HeadBlocks or manually. Their structure (count matrix, spatial coordinates and optional region annotations) is validated (`validate_input`) before downstream analyses.

2. **Region annotations**. Provided in GeoJSON format.

3. **Cell type annotations**. Automatically generated through externally annotated references (`ingest`), manually annotating clusters (default) or providing externally generated annotations.

The only input needed for SpaceBlocks is the actual ST data (1).

SpaceBlocks provides tools to ease region annotation via QuPath (ideal if collaborating with pathologists) and cell type annotation (either manual, semi-automatic or external).

## Quickstart

See the tutorial on [how to streamline a full run](docs/getting-started.md), and an [example with public data](docs/reproduction.md).

> [!TIP]
> Any run can be fully reproduced end-to-end only with the experimental design and cell type annotations.
> Color palettes for sample metadata, cell type annotations and spatial niches are fully customizable from `config/config.yaml`.

SpatialBlocks requires [Snakemake](https://snakemake.readthedocs.io) ≥ 8 and Conda/Mamba. The **Visium HD** HeadBlock
additionally needs an external [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger) ≥ 4.0.1 installation.

```bash
# 1. Get the workflow  (or: snakedeploy deploy-workflow cbib/SpaceBlocks spaceblocks --tag <version>)
git clone https://github.com/cbib/SpaceBlocks && cd SpaceBlocks

# 2. Configure
#    config/visiumhd_samples.csv.csv   — if running mode : visiumhd
#    config/core_samples.tsv   — one row per sample (+ any design columns)
#    config/config.yaml        — mode: visiumhd | xenium5k | decoupled, and paths

# 3. Check the plan
snakemake -n --sdm conda

# 4. Make the QuPath images, annotate regions in QuPath, then run the core
snakemake qupath_images       --sdm conda   # → annotate, export GeoJSON to config["geojson_path"]
snakemake run_preprocessing   --sdm conda   # validate → QC/normalise → cluster
snakemake run_postprocessing  --sdm conda   # annotate → integrate → DE → reports
```

On SLURM, add `--profile profiles/default`. Full configuration reference: [`config/README.md`](config/README.md).

## Documentation

Full, searchable documentation lives at **https://cbib.github.io/SpaceBlocks/**:

- **Design & architecture** — HeadBlock/CoreBlock split, the contract, and output tree
- **Rule reference** — every rule explained briefly, grouped by phase
- **Configuration** — explanation of all config keys and sample sheets
- **QuPath annotation** — draw and export the region annotations
- **Environments** — the Conda environments and how to lock them
- **Get started: recommendations for a full run** — [how to streamline a full run](docs/getting-started.md)
- **[Reproduce example run](docs/reproduction.md)** — example use on public data

## Citing SpaceBlocks

If you use SpaceBlocks in your research, please cite:

> _Authors. SpaceBlocks: \<title\>. \<preprint / journal\>, \<year\>. \<DOI\>_

<!-- A BibTeX entry and a Zenodo DOI badge will be added on the first release. -->

## License

Released under the MIT License — see [LICENSE](LICENSE).
