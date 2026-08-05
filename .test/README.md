# Canonical test fixture

This directory contains the single maintained integration fixture used for local
validation, GitHub Actions, and Snakemake Workflow Catalog graph construction. It
runs SpaceBlocks in `mode: decoupled` against two tiny synthetic contract h5ads,
so workflow parsing and DAG construction do not require Space Ranger, Conda
environment creation, or access to a cluster.

## Layout

```text
.test/
├── config/       Catalogue/CI configuration and input TSV files
├── data/         Deterministic synthetic S1.h5ad and S2.h5ad contracts
├── resources/    Fixture-generation utilities
└── README.md
```

All paths in `config/config.yaml` are relative to `.test`, which is the working
directory for every validation command.

## Regenerate the synthetic data

Install the generator dependencies (`anndata`, `numpy`, and `scipy`), then run
from any directory:

```bash
python .test/resources/generate_test_data.py
```

The script resolves its output directory from its own location, uses a fixed
random seed, and overwrites `.test/data/S1.h5ad` and `.test/data/S2.h5ad`.

## Validate the workflow

Run these commands from the repository root:

```bash
snakemake -s workflow/Snakefile -d .test --lint --workflow-profile none
snakemake -s workflow/Snakefile -d .test -n --workflow-profile none
snakemake -s workflow/Snakefile -d .test --rulegraph --workflow-profile none \
  | dot -Tsvg > pipeline_rulegraph.svg
```

The final command mirrors the rule-graph construction used for catalogue
presentation and requires Graphviz.

## Maintenance

- Keep `.test` as the only committed integration fixture; do not recreate
  `.tests` or duplicate its config/data elsewhere.
- Add or change sample-sheet inputs only under `.test/config`.
- Regenerate the h5ads whenever the input contract changes, and commit both the
  generator change and regenerated files.
- Keep all paths in the fixture relative to `.test` so local, CI, and catalogue
  execution remain identical.
- Run lint, dry-run, and rule-graph generation after workflow or fixture changes.
