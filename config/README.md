# Configuration

## `config.yaml` (repository root)

Main pipeline configuration.  All paths, analysis parameters, and resource
allocations are defined here.  See `workflow/schemas/config.schema.yaml` for
the full specification with types and defaults.

### Required paths

| Key | Description |
|-----|-------------|
| `samples` | Path to `config/samples.csv` |
| `snakemake_cell_markers` | Path to the cell-type marker TSV |
| `probe_set` | Visium HD probe set CSV from 10x Genomics |
| `transcriptome` | Cell Ranger / Space Ranger reference transcriptome |
| `spaceranger` | Path to the locally installed `spaceranger` executable |
| `geojson_path` | Directory with QuPath GeoJSON annotations |
| `spaceranger_processing_outdir` | Space Ranger output directory |
| `post_processing_outdir` | Post-processing output directory |
| `logdir` | Log and benchmark directory |

### Analysis parameters (`analysis:` block)

All analysis thresholds are configurable — no code edits needed.  Key
parameters include QC filters (`min_counts`, `min_cells`, `min_genes`),
feature selection (`n_top_genes`), clustering (`leiden_resolution`,
`n_neighbors`, `n_pcs`), and the resolution scan range for silhouette
evaluation.  See the schema file for full details and defaults.

### Resources (`resources:` block)

Per-rule resource allocations use SLURM-native units: `mem_mb` (integer,
megabytes), `runtime` (integer, minutes), and `threads` (integer).

## `samples.csv`

One row per Visium HD sample.  Required columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique identifier (must match FASTQ prefixes and image names) |
| `fastq_dir` | Root directory with `fastqs/` and `images/` subdirectories |
| `slide` | Visium HD slide serial number |
| `area` | Capture area (`A1`, `B1`, `C1`, or `D1`) |

Optional columns for re-sequenced libraries: `resequenced_dir`,
`new_resequenced_dir`, `third_resequencing_dir`.

## `snakemake_cell_markers.tsv`

A column-oriented TSV where each column header is a cell-type label and the
rows below it list marker gene symbols.  Columns can have different lengths
(shorter columns are padded with empty strings).
