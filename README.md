# Visium HD Snakemake Pipeline

End-to-end [Snakemake](https://snakemake.github.io/) workflow for 10x Genomics
**Visium HD** spatial transcriptomics — from raw FASTQ files to cell-level
clustering, spatial annotation, and cell-type marker visualisation.

## Pipeline overview

```
FASTQ + images ─► Space Ranger count ─┬─► generate_qupath_tiff
                                       │       ↓
                                       │   (user annotates in QuPath,
                                       │    exports GeoJSON)
                                       │       ↓
                                       └─► post_processing (bin2cell)
                                               ├─ QC filtering
                                               ├─ QuPath region annotations
                                               ├─ Seurat normalisation + Harmony
                                               ├─ UMAP + Leiden clustering
                                               ├─ Silhouette + clustree evaluation
                                               ├─ Differential expression
                                               ├─ Cell-type marker visualisation
                                               └─ AnnData export (.h5ad)
```

## Repository layout

```
├── .gitignore
├── .snakemake-workflow-catalog.yml     # Workflow Catalog metadata
├── config.yaml                         # Main pipeline configuration
├── config/
│   ├── README.md                       # Configuration documentation
│   ├── samples.csv                     # Sample sheet
│   └── snakemake_cell_markers.tsv      # Cell-type marker genes
├── profiles/
│   └── default/
│       └── config.yaml                 # SLURM execution profile
├── .tests/                             # Minimal test fixtures
│   ├── config.yaml
│   └── config/
│       ├── samples.csv
│       └── snakemake_cell_markers.tsv
└── workflow/
    ├── Snakefile                       # Entry point (with schema validation)
    ├── envs/
    │   └── bin2cell.yaml               # Conda environment (minimal)
    ├── rules/
    │   ├── spaceranger.smk             # Space Ranger count
    │   ├── generate_qupath_tiff.smk    # Low-res TIFF for QuPath annotation
    │   └── post_processing.smk         # bin2cell analysis
    ├── schemas/
    │   ├── config.schema.yaml          # JSON Schema for config.yaml
    │   └── samples.schema.yaml         # JSON Schema for samples.csv
    └── scripts/
        ├── generate_qupath_tiff.py     # TIFF generation script
        └── post_processing.py          # Main analysis script
```

## Requirements

- **Snakemake ≥ 8.0** with `snakemake-executor-plugin-slurm`
- **Space Ranger 4.0.1** installed locally (licence does not permit
  redistribution in containers; download from
  [10x Genomics](https://www.10xgenomics.com/support/software/space-ranger/downloads))
- **Conda / Mamba** (for automatic environment creation)

## Quick start

```bash
# 1. Edit config.yaml — set spaceranger path and all directory paths
vim config.yaml

# 2. Edit config/samples.csv with your samples
vim config/samples.csv

# 3. Dry run
snakemake --profile profiles/default -n

# 4. Run Space Ranger + generate QuPath TIFFs
snakemake --profile profiles/default -- generate_all_qupath_tiffs

# 5. Annotate TIFFs in QuPath, export GeoJSON to geojson_path directory

# 6. Run full pipeline (post-processing will pick up GeoJSON annotations)
snakemake --profile profiles/default
```

## QuPath annotation workflow

The pipeline includes a dedicated rule (`generate_qupath_tiff`) that creates
low-resolution TIFF images suitable for annotation in QuPath:

1. **Run the rule** — after Space Ranger completes, invoke:
   ```bash
   snakemake --profile profiles/default -- generate_all_qupath_tiffs
   ```
   This creates `he_{sample_id}_mpp{mpp}.tiff` in the `geojson_path` directory.

2. **Annotate in QuPath** — open the TIFF, draw region annotations
   (e.g. tumour, stroma, necrosis).

3. **Export as GeoJSON** — save the annotation file as
   `he_{sample_id}_mpp{mpp}.geojson` in the same `geojson_path` directory.

4. **Run post-processing** — the pipeline will automatically detect and
   integrate the GeoJSON annotations.  If no GeoJSON is found for a sample,
   the analysis continues without region labels.

## Configuration

All configuration lives in `config.yaml`.  See `config/README.md` for full
documentation, and `workflow/schemas/` for the validation schemas.

Key design decisions:

- **All analysis thresholds are in `config.yaml`** — QC filters, HVG count,
  kNN parameters, Leiden resolution, silhouette scan range, DE gene count.
  No code edits needed to change parameters.

- **Schema validation** — both `config.yaml` and `samples.csv` are validated
  at DAG-construction time, catching typos before any jobs run.

- **Real output tracking** — Space Ranger outputs (`web_summary.html`,
  `raw_feature_bc_matrix.h5`) and post-processing deliverables
  (`adata_SH_*.h5ad`, `bdata_SH_*.h5ad`, `cluster_markers_*.tsv`) are
  tracked as actual Snakemake outputs.

- **Minimal conda environment** — the env YAML specifies only the direct
  dependencies at major.minor level, letting the solver resolve the full
  tree.  This is smaller, more portable, and faster to build than a pinned
  export of 240+ packages.

## Outputs

Per sample, the pipeline produces:

| Output | Location |
|--------|----------|
| Space Ranger web summary | `{sr_outdir}/{sample}/outs/web_summary.html` |
| QuPath TIFF | `{geojson_path}/he_{sample}_mpp{mpp}.tiff` |
| Processed AnnData | `{pp_outdir}/{sample}/adata_SH_{sample}.h5ad` |
| HVG-subset AnnData | `{pp_outdir}/{sample}/bdata_SH_{sample}.h5ad` |
| Cluster markers TSV | `{pp_outdir}/{sample}/cluster_markers_{sample}.tsv` |
| QC plots | `{pp_outdir}/{sample}/QC_plots/` |
| UMAP / spatial plots | `{pp_outdir}/{sample}/` |
| Clustering evaluation | `{pp_outdir}/{sample}/clustering_evaluation/` |
| Cell-type marker plots | `{pp_outdir}/{sample}/cell_type_markers/` |

## Testing

```bash
snakemake --profile profiles/default -n \
    --configfile .tests/config.yaml \
    --config samples=.tests/config/samples.csv \
           snakemake_cell_markers=.tests/config/snakemake_cell_markers.tsv
```

## Linting

```bash
snakemake --lint -s workflow/Snakefile
```
