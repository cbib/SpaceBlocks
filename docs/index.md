# SpaceBlocks

**A Snakemake workflow for single-cell-resolution spatial transcriptomics.**

SpaceBlocks turns raw spatial data from multiple platforms into one standardized, analysis-ready
object (the *contract* h5ad) and runs a full, reproducible analysis on top of it — quality control,
clustering, annotation, integration, differential expression, and spatial neighbourhood/niche
analysis — identically on a workstation or an HPC cluster.

<p align="center"><img src="img/overview.svg" alt="SpaceBlocks workflow overview" width="760"></p>

## How it fits together

A **Headblock** is technology-specific and builds the standardized contract h5ad; the
**Coreblock** is technology-agnostic and consumes it. Adding a new platform means writing a new
head that emits the contract — the analysis core never changes.

- **[Pipeline & architecture](pipeline.md)** — the Headblock/Coreblock split, the contract spec, the output tree.
- **[Rule reference](rules.md)** — every rule, grouped by phase.
- **[Configuration](configuration.md)** — modes, sample sheets, and all config keys.
- **[Environments](environments.md)** — the per-rule Conda environments and how to lock them.
- **[Reproduce the paper run](reproduction.md)** — full validation on public data.

## Quickstart

See the [repository README](https://github.com/<owner>/SpaceBlocks#quickstart) for the minimal
run, or the [Configuration](configuration.md) page for the full setup.
