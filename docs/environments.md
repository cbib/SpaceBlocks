# Index

- [Conda environments & version pinning](#conda-environments--version-pinning)
  - [Generating exact pins from the Snakemake-managed environments](#generating-exact-pins-from-the-snakemake-managed-environments)
  - [Alternative: conda-lock (multi-platform locks)](#alternative-conda-lock-multi-platform-locks)
  - [Note on the DE environment](#note-on-the-de-environment)

# Conda environments & version pinning

The pipeline provisions one conda environment per rule group (run with `--use-conda`):

| env file | used by | pinning |
|---|---|---|
| `envs/visiumhd.yaml` | Headblock + most Coreblock rules (scanpy/squidpy stack) | `>=` lower bounds |
| `envs/pseudobulk_aggregate.yaml` | pseudobulk aggregation, gene exploration (decoupler) | `>=` lower bounds |
| `envs/pseudobulk_de.yaml` | pseudobulk DE (R / DESeq2) | pinned to Bioconductor 3.18 / R 4.3 |
| `envs/spatial_niches.yaml` | `spatial_niches` rule | `>=` lower bounds |
| `envs/xenium5k.yaml` | xenium5k HeadBlocks | `>=` lower bounds |


`>=` bounds are reproducible enough for day-to-day use but not for archival reproducibility (a future solve may pick newer, potentially breaking versions). For a publication release, we recommend to generate and share the **exact** locks from the environments you actually tested, as below.

> [!TIP]
> For a reproducible release, commit one **lock file per environment** next to the
> `envs/*.yaml`. Generate them from the environments you actually tested (see below),
> e.g. `conda list -p .snakemake/conda/<hash>_ --explicit > envs/<env>_linux-64.lock`,
> and use `conda create --file envs/<env>_linux-64.lock` to rebuild the exact stack.
> Keep the loose `envs/*.yaml` as the maintainable/reusable surface and the locks as
> the exact-reproducibility surface.

## Generating exact pins from the Snakemake-managed environments

Snakemake builds each env under `.snakemake/conda/<hash>_/` on first `--use-conda` run.

```bash
# 1. Create the envs without running the workflow (or just run it once).
snakemake --use-conda --conda-create-envs-only --cores 1

# 2. List each env file and the prefix Snakemake built for it.
snakemake --use-conda --list-conda-envs
#   ... prints:  workflow/envs/<envname>.yaml   .snakemake/conda/<hash>_

# 3. Export EXACT versions from a given env (choose one style):

#  a) versions only, cross-platform (recommended for envs/*.pinned.yaml):
conda env export -p .snakemake/conda/<hash>_ --no-builds \
  > workflow/envs/<envname>.pinned.yaml

#  b) fully reproducible, platform-specific (versions + builds + URLs):
conda list -p .snakemake/conda/<hash>_ --explicit \
  > workflow/envs/<envname>.linux-64.lock
#     recreate with:  conda create --name X --file workflow/envs/<envname>.linux-64.lock
```

Repeat per env, then point the rules at the `.pinned.yaml` files (or keep the loose
`.yaml` for development and the `.lock`/`.pinned.yaml` for releases).

## Alternative: conda-lock (multi-platform locks)

`conda-lock` resolves once and pins every transitive dependency with hashes, giving the strongest reproducibility guarantee across machines.

```bash
pip install conda-lock
conda-lock lock -f workflow/envs/<envname>.yaml -p linux-64 -p osx-64 --kind explicit
#   -> conda-linux-64.lock, conda-osx-64.lock
```

## Note on the DE environment

`workflow/envs/pseudobulk_de.yaml` is pinned to a coherent Bioconductor 3.18 / R 4.3 release. If your working environment used a different Bioconductor/R release, replace it with an exact export (above) so the published env matches what you validated.

Notice that the minor-version pins here are a sensible default, not a guarantee that they equal your tested versions.
