# Index

- [Outputs](#outputs)
    - [Output directory structure](#output-directory-structure)
- [Preprocessing outputs](#preprocessing-outputs)
- [Postprocessing outputs](#postprocessing-outputsç)
- [Exploration outputs](#exploration-outputs)

# Outputs

What each step of the workflow writes, and why it is useful. Everything lands under
`post_processing_outdir` (set in `config.yaml`); per-sample results live in `Samples/<sample>/`.

## Output directory structure

The tag after each folder marks which CoreBlock step produces it — **(pre)** preprocessing,
**(post)** postprocessing, **(expl)** exploration.

```
<post_processing_outdir>/
├── Samples/<sample>/
│   ├── <sample>_unfiltered.h5ad                    (contract — head or decoupled input)
│   ├── validation/                                 (pre)   input_validation.json
│   ├── qc_sweep/                                    (pre)   optional pre-filter QC diagnostics
│   ├── adata_<sample>.h5ad, metadata_<sample>.tsv   (pre)
│   ├── <sample>_report.tsv                          (pre)
│   ├── leiden_resolution_<r>/                        (pre)
│   ├── adata_<sample>_ingested.h5ad, ingest/         (pre)   optional
│   ├── adata_<sample>_annotated.h5ad, annotation/    (post)
│   └── neighbourhood_analysis/<annot_type>/          (post)
├── spatial_niches/                                  (pre)   optional
├── cluster_annotations_template.tsv                 (pre)   fill this in, then annotate
├── integrated_samples/                              (post)  concatenated + integrated objects, report
├── pseudobulk/<annot_type>/<analysis_level>/        (post)  aggregated/ + de_results/
├── Subcompartments/<name>/                          (post)  optional
└── gene_exploration/                                (expl)
```

## Preprocessing outputs

| Rule | File | What it is |
| --- | --- | --- |
| `validate_input` | `validation/input_validation.json` | Confirms the contract is valid and records what was found |
| `qc_sweep` | `qc_sweep/qc_*.png` | Shows where candidate QC thresholds would remove cells, before filtering |
| `preprocess_umap` | `adata_<sample>.h5ad` | Filtered, normalised object with PCA, UMAP and multi-resolution clusters |
| `preprocess_umap` | `metadata_<sample>.tsv` | Per-cell table of clusters and QC metrics |
| `preprocess_umap` | `<sample>_report.tsv` | Cells before vs after filtering and the cut-offs applied |
| `leiden_analysis` | `leiden_resolution_<r>/` | UMAP and spatial plots for each clustering resolution |
| `generate_annotation_template` | `cluster_annotations_template.tsv` | Empty cluster to cell-type sheet for you to fill in |
| `ingest_ref` | `adata_<sample>_ingested.h5ad`, `ingest/` | Optional reference-based cell-type labels and plots |
| `spatial_niches` | `spatial_niches/tsv/niche_<sample>.tsv`, `plots/` | Optional spatial niches per cell and their plots |

## Postprocessing outputs

| Rule | File | What it is |
| --- | --- | --- |
| `annotate_cells` | `adata_<sample>_annotated.h5ad` | Object with final cell-type labels applied |
| `annotate_cells` | `annotation/` | Per-sample annotation and marker plots |
| `integrate_samples` | `integrated_samples/harmony_integrated.h5ad` | All samples concatenated and batch-integrated |
| `integrate_samples` | `integrated_samples/sketched.h5ad` | Geometric-sketch subsample for fast plotting |
| `pseudobulk_aggregate` | `pseudobulk/<annot_type>/<analysis_level>/aggregated/` | Pseudobulk count matrices and their PCA QC |
| `pseudobulk_de` | `pseudobulk/<annot_type>/<analysis_level>/de_results/` | Differential-expression tables and plots |
| `neighbourhood_analysis` | `neighbourhood_analysis/<annot_type>/` | Cell-type co-occurrence and enrichment per sample |
| `subcluster` | `Subcompartments/<name>/` | Optional re-clustering of a chosen compartment |
| `sample_report` | `integrated_samples/samples_report.pdf` | One page per sample with UMAP, spatial, dotplot and barplot |

## Exploration outputs

| Rule | File | What it is |
| --- | --- | --- |
| `explore_genes_integrated` | `gene_exploration/expression_ranges.tsv` | Shared expression ranges so plots use one colour scale |
| `explore_genes_integrated` | `gene_exploration/<entry>/Integrated/` | Signature scores and dotplots on the integrated object |
| `explore_genes_sample` | `gene_exploration/<entry>/Spatial/` | Spatial maps of each gene or signature per sample |
| `explore_genes_sample` | `gene_exploration/<entry>/Dotplots/` | Per-sample dotplot composites |
