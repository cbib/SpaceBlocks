# Reproduce the paper run

The small case in `.tests/` exists only to exercise the DAG in CI (and to render the workflow's
tube map on the Snakemake catalog). It is **not** a scientific validation.

The full, public-data validation used for the publication lives under **`reproduction/`** in the
repository. The datasets themselves are **not** committed — they are downloaded from their public
source at run time — so the repository stays small and catalog-clonable.

```bash
cd reproduction
bash fetch_data.sh                       # download the public datasets (GEO / Zenodo / 10x)
snakemake --sdm conda \
  --configfile config/config.yaml \
  --directory .                          # run the full workflow on the fetched data
```

| Path | What it is |
| --- | --- |
| `reproduction/config/` | full-run `config.yaml` + `core_samples.tsv` pointing at the public data |
| `reproduction/fetch_data.sh` | downloads the datasets into a git-ignored `reproduction/data/` |
| `reproduction/run.sh` | the exact commands used for the published run |

See [Configuration](configuration.md) for what each key means.
