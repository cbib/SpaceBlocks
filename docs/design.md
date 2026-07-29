# Index

- [Design & architecture](#design--architecture)
- [Repository structure](#repository-structure)
- [The Blocks](#the-blocks)
- [The validation contract](#the-validation-object-contract)
- [Output structure](#output-structure)
- [Key design decisions](#key-design-decisions)
- [Extending SpaceBlocks](#extending-spaceblocks)
  - [Extending HeadBlocks to a new platform](#extending-headblocks-to-a-new-platform)
  - [Extending the CoreBlocks](#extending-the-coreblocks)

# Design & architecture

SpaceBlocks is built around one idea: a technology-specific **HeadBlock** turns raw data into a standardized *contract* h5ad AnnData object, and a technology-agnostic **CoreBlock** does all the analysis on top of it.

Adding a platform means writing a new HeadBlock; the CoreBlock is shared and stable across technologies.

Importantly, the CoreBlocks immediately admit external data, as soon as they pass the `validate_input` rule (see [below](design.md#the-validation-contract)).

## Repository structure

```
SpaceBlocks/
├── workflow/                 the executable workflow (auto-found by Snakemake)
│   ├── Snakefile             globals, mode selection, targets
│   ├── rules/*.smk           one file per rule (+ common.smk for shared helpers)
│   ├── scripts/*.py, *.R     rule implementations
│   ├── envs/*.yaml           one Conda env per rule group
│   └── schemas/*.yaml        config + sample-sheet validation
├── config/                   config.yaml, README.md (config reference), sample sheets
├── docs/                     this documentation site (MkDocs)
├── .tests/                   tiny synthetic decoupled dataset for CI
├── reproduction/             full public-data runs (data fetched, not committed)
├── tools/                    stand-alone helper scripts
└── profiles/default/         SLURM profile (retries, resources)
```

## The Blocks

```
── HeadBlock (technology-specific) ──────────────────────────────────────────────
  mode: visiumhd      spaceranger_count_vhd → generate_qupath_vhd → prepare_input_vhd
  mode: xenium5k      convert_zarr_x5k      → generate_qupath_x5k  → prepare_input_x5k
  mode: decoupled     (no HeadBlock — you provide the contract h5ads)
        │  ‖ COREBLOCK INPUT POINT (validate_input): contract h5ad structure validation
        ▼  each HeadBlock writes the unfiltered contract h5ad; although these can also be prepared manually (mode: decoupled)
── CoreBlock (technology-agnostic) ──────────────────────────────────────────────
  STEP 1 · preprocessing
    validate_input → [qc_sweep] → preprocess_umap → leiden_analysis
                                 → generate_annotation_template
                                 → ingest_ref (opt.) → spatial_niches (opt.)
  STEP 2 · postprocessing
    annotate_cells → integrate_samples → pseudobulk_aggregate → pseudobulk_de
                   → neighbourhood_analysis → subcluster → sample_report
  STEP 3 · exploration
    explore_genes_integrated → explore_genes_sample
```

See the [Rule reference](rules.md) for a one-paragraph description of every rule.

## The validation object contract

The point where the HeadBlocks meet the CoreBlocks is the `validate_input` rule. This rule takes as an **input** the standardized h5ad AnnData files (**the *contracts***, generated manually or via the HeadBlocks) and **outputs** a `.json` validation file, which confirms the data is suitable for entering pre- and postprocessing steps.

It is a single coupling point with one h5ad per sample. Every HeadBlock produces the same contract file structure, so the CoreBlock is identical for all technologies:

- **Path** — one config key (`contract.unfiltered_h5ad`, derived in the Snakefile), *written* by
  the active HeadBlock's `prepare_input_*` and *read* by `validate_input`, `qc_sweep`, and
  `preprocess_umap`, so producer and consumers can never drift.
- **Shape of AnnData object:**
  - `X` — raw integer counts (no filtering, no normalisation).
  - `obsm["spatial"]` — per-cell `(x, y)` centroids, finite.
  - `uns["spatial"][sample]` — `{images: {hires}, scalefactors}` when a tissue/morphology image
    exists; omitted otherwise (the CoreBlock then scatters on coordinates).
  - `obs` — `sample` (required); `region_annotation` when a GeoJSON was provided; and
    `cell_id` (optional but strongly recommended, and **must equal `obs_names`** to prevent
    unexpected behaviour, since downstream joins key on it — `validate_input` warns otherwise).
- **The image is a HeadBlock concern.** The CoreBlock never rebuilds an image — it consumes the
  embedded `uns["spatial"]` image when present and otherwise scatters on `obsm["spatial"]`.

For a worked walk-through of building these files by hand (naming, layout, config), see
[Preparing inputs for decoupled mode](reproduction.md#preparing-inputs-for-decoupled-mode).

## Output structure

```
{post_processing_outdir}/
├── Samples/{sample}/
│   ├── {sample}_unfiltered.h5ad                 CONTRACT H5AD (from HadBlocks: prepare_input_{vhd,x5k})
│   ├── QuPath_image/…                            generate_qupath_{vhd,x5k}
│   ├── validation/input_validation.json         validate_input
│   ├── qc_sweep/                                 qc_sweep (optional)
│   ├── adata_{sample}.h5ad, metadata_{sample}.tsv, {sample}_report.tsv   preprocess_umap
│   ├── leiden_resolution_*/                      leiden_analysis
│   ├── annotation/, adata_{sample}_annotated.h5ad   annotate_cells
│   └── neighbourhood_analysis/{annot_type}/      neighbourhood_analysis
├── integrated_samples/                           concatenated / harmony / sketched / samples_report.pdf
├── pseudobulk/{annot_type}/{analysis_level}/     aggregated/ + de_results/
├── Subcompartments/{name}/{Harmony,NoHarmony}/   subcluster
├── spatial_niches/                               tsv/, plots/
├── gene_exploration/                             expression_ranges.tsv, {entry}/…
└── cluster_annotations_template.tsv              generate_annotation_template
```

## Key design decisions

- **HeadBlock/CoreBlock split.** A technology-specific set of rules (HeadBlock) produces a standardized unfiltered AnnData file (contract); the common CoreBlocks consume it for the analyses. Inclusion of new platform/s only requires the development of a new HeadBlock, and the CoreBlock does not need to be changed, unless it has to be expanded.
- **`validate_input` is a DAG gate.** The division in ensured by a validation rule, which can pass or give a hard/soft failure.
  - A `.json` is written if the contract structure is validated.
  - A hard failure prevents the `.json` from being written, and thus the CoreBlock from running.
  - Soft issues (missing region annotation, no mito genes, no image) are recorded, but do not necessarily prevent CoreBlocks from running (config file has a paramter for soft-passing).
- **Coherent naming.** Head rules carry a 3-letter technology code (`_vhd`, `_x5k`) so the organization is easy to follow and heads can coexist.
- **`qc_sweep` rule is diagnostic only** — it never filters, clusters, or writes an h5ad.
- **External annotation takes over.** When enabled, it becomes the primary annotation everywhere. The [Configuration](configuration.md) allows flexibility to retain all cells or remove externally unannotated ones.
- **Config-driven colours** — regions, sample metadata, and cell types, applied consistently across every plot, with a grey fallback for undefined levels. This allows precise and consitent color representations through the analyses.
- **Retries scale memory.** `mem_mb` grows with the attempt number, so an OOM-killed job is resubmitted with more RAM.

## Extending SpaceBlocks

The aim of SpaceBlocks is to allow for long-term maintainable and extensible Spatial Transcriptomics analyses across platforms and tools. Below you can find guidelines for development.

If you want to contribute or extend SpaceBlocks, see [CONTRIBUTING.md](https://github.com/cbib/SpaceBlocks/blob/main/CONTRIBUTING.md) for the conventions and validation workflow.

### Extending HeadBlocks to a new platform

If you wished to run SpaceBlocks but your platform is not listed among the available SpaceBlocks HeadBlocks, there are two options available:
1. Running SpaceBlocks in `mode: decoupled` (**easier, recommended**).
2. Writting a new SpaceBlocks HeadBlock.

See [uncoupled demos](reproduction.md) for examples about how to run SpaceBlocks in decoupled mode.

If you wished to write a new HeadBlock:

- Create a 3-letter `<code>` to identify the technology.
- Write the Snakemake rules and Python scripts to prepare the AnnData h5ad contracts. Name them as `rules/*_<code>.smk` and `scripts/*_<code>.py`.
- Add the mode in the `Snakefile` parameter `KNOWN_HEADS`.
- Create a `demo/` (public data download, GeoJSON region annotations, short documentation) with a prepared config example and instructions for reproducibility with a test run.
- Include the documentation in the `reproduction.md`.
- Create a Pull Request to the `development` branch. We will revise it promptly.

### Extending the CoreBlocks

SpaceBlocks includes many utilities for ST analyses, but there is room for extension and inclusion of flexible analysis parameters (i.e. methods for normalization, integration, variable gene selection, etc.).

If you wished to contribute, we will be happy to revise contributions to the SpaceBlocks workflow.
