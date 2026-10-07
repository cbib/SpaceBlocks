<p align="center">
  <img src="docs/img/spaceblocks_logo.svg" alt="SpaceBlocks" width="570">
</p>
<p align="center">
  <a href="https://cbib.github.io/SpaceBlocks/"><img src="https://img.shields.io/badge/docs-online-blue" alt="Documentation"></a>
  <a href="https://github.com/cbib/SpaceBlocks/actions/workflows/tests.yml"><img src="https://github.com/cbib/SpaceBlocks/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://github.com/cbib/SpaceBlocks/blob/main/LICENSE.md"><img src="https://img.shields.io/github/license/cbib/SpaceBlocks" alt="License"></a>
  <a href="https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/cbib/SpaceBlocks.html"><img src="https://img.shields.io/badge/Snakemake-workflow-039475" alt="Snakemake workflow"></a>
</p>

<p align="center">
  <a href="https://cbib.github.io/SpaceBlocks/">Read the online documentation</a> ·
  <a href="https://cbib.github.io/SpaceBlocks/demos/">Demos</a> ·
  <a href="CONTRIBUTING.md">Contributing</a>
</p>

## Streamlined, reproducible, spatial transcriptomics analyses with Snakemake

SpaceBlocks is a technology-agnostic Snakemake workflow for **expert-supervised cohort-level analysis of single cell resolution spatial transcriptomics (ST)**.

It supports scalable and reproducible quality control, clustering, annotation, integration, differential expression, and spatial neighbourhood and niche analysis.

Its modular architecture supports both public and in-house ST datasets while remaining maintainable and extensible.

## Quickstart glossary

- **HeadBlocks** — technology-specific sets of rules that prepare the raw data files into a contract h5ad.
- **CoreBlocks** — technology-agnostic sets of rules that streamline ST processing steps.
- **Contract** — The standardized AnnData h5ad files, validated before running the CoreBlocks.


## Overview

The workflow is divided into (optional) technology-specific **HeadBlocks** and common **CoreBlocks** that streamline preprocessing, postprocessing, and informative exploration of results.

<br/>

<p align="center">
  <img src="images/main.png"
       alt="SpaceBlocks workflow overview"
       width="700">
</p>

> **Overview of the SpaceBlocks workflow**. **(a)** Snakevision tubemap scheme showcasing the workflow rules, coloured by technology and processing step. The AnnData contract validation entry point and user-supervised steps are displayed in **bold font**. **(b)** Schematic representation of the outputs produced by the CoreBlocks subworkflow and main design features of the pipeline. **(c)** Runtime benchmarking across Visium HD, Xenium 5K and MERSCOPE datasets. β denotes the fitted scaling exponent; doubling the cell count multiplies the estimated wall-clock time by 2^β. **(d)** SpaceBlocks processing of the multiple sclerosis MERSCOPE dataset from Feng et al. (2025). Left, BANKSY niches annotated in the publication, shown on a UMAP. Right, the number of unique domain-specific upregulated genes identified by one-vs-rest Wald tests (Benjamini–Hochberg adjusted p-value < 0.05 and |log2FC| > 0.25).

<br/>

SpaceBlocks consumes ST data and can additionally use region and cell-type annotations:

1. **ST formatted AnnData objects**. Either generated from the HeadBlocks, or manually formatted as a standardized h5ad AnnData object. For brevity, **we refer to each of these objects as THE CONTRACT**. Their structure (count matrix, spatial coordinates and optional region annotations) is validated (`validate_input`) before downstream analyses.

2. **Region annotations (optional)**. Provided as GeoJSON files in HeadBlock modes, or already
   stored in `obs["region_annotation"]` in a decoupled contract.

3. **Cell type annotations**. Automatically transferred from a reference (`ingest`),
   manually assigned from clusters (default), or supplied as external labels.

SpaceBlocks provides tools to ease region annotation via QuPath (ideal if collaborating with pathologists) and cell type annotation (either manual, ingested or external).

## Quickstart

SpaceBlocks requires [Snakemake](https://snakemake.readthedocs.io) and Conda/Mamba. The pipeline was developed under Snakemake v9.13.7. The **Visium HD** HeadBlock additionally needs an external [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger) ≥ 4.0.1 installation.

If you are not familiar with Snakemake, there are three key files to configure:

1. `workflow/Snakefile` — loads the config file/s (normally, at `config/config.yaml`), and has the instructions to execute the workflow. It is normally static, so you should not modify it.
2. The `config.yaml` — defines the parameters for your workflow run (i.e. where the input and output should be found, whether to use an external reference for cell type annotation, etc.). Copy `config/config.yaml.template` to `config/config.yaml`, then configure it for your project. The config may require other files needed for the run (in the case of SpaceBlocks, `config/core_samples.tsv`).
3. The `profiles/<profile>/config.yaml` — Snakemake is a workflow manager that allows parallelization in systems with job schedulers (i.e. `Slurm`), or local execution. **The profile needs to be adapted to your system.** We provide a default profile that runs the `test` pipeline locally and a slurm profile we used on our HPC.

The full [configuration reference](https://cbib.github.io/SpaceBlocks/configuration/) presents a detailed explanation on how to set `config.yaml`.

We recommend reading the [configuration reference](https://cbib.github.io/SpaceBlocks/configuration/) first, and then starting with the tutorial on [how to streamline a full run](https://cbib.github.io/SpaceBlocks/getting-started), and an [example with public data](https://cbib.github.io/SpaceBlocks/demos/).

> [!TIP]
> Once the input data is available, a run can be reproduced from its configuration, sample sheet, and any manual or external annotations.
>
> Color palettes for sample metadata, cell type annotations and spatial niches are fully customizable from `config/config.yaml`.

In short, to run the workflow:

```bash
# 1. Get the workflow  (or: snakedeploy deploy-workflow cbib/SpaceBlocks spaceblocks --tag <version>)
git clone https://github.com/cbib/SpaceBlocks && cd SpaceBlocks

# 2. Configure the config/config.yaml and the Snakemake profile in profiles/
cp config/config.yaml.template config/config.yaml
#    config/visiumhd_samples.csv   — if running mode : visiumhd
#    config/core_samples.tsv       — one row per sample (+ any design columns)
#    config/config.yaml            — mode: visiumhd | xenium5k | atera | merscope | decoupled
#    profiles/**/config.yaml       — local and Slurm execution examples
#
#    Note: Invoke profile manually (e.g. --profile profiles/default) or
#          set environment variable (e.g. export SNAKEMAKE_PROFILE=profiles/default )

# 3. Dry-run your SpaceBlocks configuration
snakemake -n --sdm conda

# 4. Make the QuPath images (only when using the HeadBlocks)
snakemake qupath_images       --sdm conda
# 5. Annotate regions in QuPath and screen the QCs
snakemake  qc_sweep_all       --sdm conda    # sweep QC thresholds
# 6. Commit the QCs and run the CoreBlocks
snakemake run_preprocessing   --sdm conda    # validate → QC/normalise → cluster
snakemake run_postprocessing  --sdm conda    # annotate → integrate → DE → reports + subclusters
snakemake run_exploration     --sdm conda    # explore genes (integrated and per sample)
```

## Documentation

Full, searchable documentation lives at **https://cbib.github.io/SpaceBlocks/**:

- **[Configuration](https://cbib.github.io/SpaceBlocks/configuration/)** — explanation of all config keys and sample sheets.
- **[Get started: recommendations for a full run](https://cbib.github.io/SpaceBlocks/getting-started/)** — how to streamline a full run.
- **[Demos](https://cbib.github.io/SpaceBlocks/demos/)** — public data end-to-end example runs.
- **[QuPath annotation](https://cbib.github.io/SpaceBlocks/qupath-tutorial/)** — tutorial to easily draw and export the region annotations.
- **[Workflow design, architectural decisions and output structure](https://cbib.github.io/SpaceBlocks/design/)** — a brief explanation on the HeadBlock/CoreBlock split, the standardized *contract* h5ad AnnData object, and output tree.
- **[Rule reference](https://cbib.github.io/SpaceBlocks/rules/)** — every rule explained briefly, grouped by HeadBlock/CoreBlock phase.
- **[Outputs](https://cbib.github.io/SpaceBlocks/outputs/)** — what each step produces and why it is useful.
- **[Environments](https://cbib.github.io/SpaceBlocks/environments/)** — the Conda environments and how to lock them.

## Citing SpaceBlocks

If you use SpaceBlocks in your research, please cite:

> _Authors. SpaceBlocks: \<title\>. \<preprint / journal\>, \<year\>. \<DOI\>_

<!-- A BibTeX entry and a Zenodo DOI badge will be added on the first release. -->

## Contributing

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) for the architecture, conventions, and how to validate a change.

## License

Released under the MIT License — see [LICENSE](LICENSE.md).
