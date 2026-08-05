# SpaceBlocks

**A Snakemake workflow for single-cell-resolution spatial transcriptomics.**

This site documents SpaceBlocks, a technology-agnostic Snakemake workflow to **semi-automatically analyse single-cell-resolution spatial transcriptomics (ST).**

**New here?** Start taking a look at the **[configuration](configuration.md)**, follow with **[get started](getting-started.md)** for a complete run and check the [examples with public data](demos.md).

## Quickstart glossary

- **HeadBlocks** — technology-specific sets of rules that prepare the raw data files into a contract h5ad.
- **CoreBlocks** — technology-agnostic sets of rules that streamline ST processing steps.
- **Contract** — the standardized AnnData h5ad files, validated before running the CoreBlocks.

Implemented HeadBlocks: **Visium HD**, **Xenium 5K**, **MERSCOPE**, and **Atera** (alpha — the platform has not shipped yet; see the [rule reference](rules.md)).

## Full documentation

- **[Configuration](configuration.md)** — explanation of all config keys and sample sheets.
- **[Get started: recommendations for a full run](getting-started.md)** — how to streamline a full run.
- **[Demos](demos.md)** — public data end-to-end example runs.
- **[QuPath annotation](qupath-tutorial.md)** — tutorial to easily draw and export the region annotations.
- **[Workflow design, architectural decisions and output structure](design.md)** — a brief explanation on the HeadBlock/CoreBlock split, the contract, and output tree.
- **[Rule reference](rules.md)** — every rule explained briefly, grouped by HeadBlock/CoreBlock phase.
- **[Outputs](outputs.md)** — what each step produces and why it is useful.
- **[Environments](environments.md)** — the Conda environments and how to lock them.

---

- **Source code:** [github.com/cbib/SpaceBlocks](https://github.com/cbib/SpaceBlocks)
- **Snakemake Workflow Catalog:** [entry for SpaceBlocks](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbib/SpaceBlocks)

---
