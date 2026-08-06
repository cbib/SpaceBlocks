# Contributing to SpaceBlocks

Thanks for your interest in improving SpaceBlocks!

This guide covers the architecture you need to know, the conventions the codebase relies on, and how to validate a change before opening a PR.
For the full design rationale, see the [documentation](https://cbib.github.io/SpaceBlocks/).

## Architecture in one minute

SpaceBlocks is split into two halves that meet at one object (the full picture is in [Design & architecture](https://cbib.github.io/SpaceBlocks/design/)):

- a technology-specific **HeadBlock** (`visiumhd`, `xenium5k`, `atera`, or `merscope`) turns raw data into a **standardized AnnData object** saved as an h5ad file, which we refer to as a **contract h5ad** for brevity;
- a technology-agnostic **CoreBlock** validates that contract and runs all the analysis.

`mode: decoupled` skips the head entirely and consumes contract h5ads you provide. **The contract is the only coupling point**: get it right and everything downstream just works.
Its exact shape (raw counts in `X`, `obsm["spatial"]`, optional `uns["spatial"]` image, `obs` keys) is specified in [Design → the validation object contract](https://cbib.github.io/SpaceBlocks/design/#the-validation-object-contract), and `validate_input` enforces it before any core rule runs.

The practical consequence: **to support a new platform you write a new HeadBlock that emits the contract, and the CoreBlock is untouched.**

## Repository layout

```
workflow/
├── Snakefile         globals, mode selection, named targets
├── rules/*.smk       one file per rule (+ common.smk for shared helpers)
├── scripts/*.py,*.R  rule implementations
├── envs/*.yaml       one Conda env per rule group (+ *_linux-64.lock)
└── schemas/*.yaml    config + sample-sheet validation
config/               config.yaml, README.md (config reference), sample sheets
docs/                 the MkDocs site
.test/               tiny synthetic decoupled fixture for CI
reproduction/         public-data worked examples
```

## Development setup

```bash
git clone https://github.com/cbib/SpaceBlocks && cd SpaceBlocks
# Snakemake >= 8 and Conda/Mamba are the only host requirements.
snakemake -n --sdm conda      # dry-run: builds the DAG, validates the config, provisions envs
snakemake -s workflow/Snakefile -d .test -n --workflow-profile none  # decoupled smoke test
```

Each rule group has its own Conda env, provisioned by `--sdm conda`. Keep the loose `envs/*.yaml`
as the maintainable surface; the `*_linux-64.lock` files are the exact-reproducibility surface (see
[Environments](https://cbib.github.io/SpaceBlocks/environments/)).

## Conventions (please follow these)

These are load-bearing, and most past bugs we experienced during development came from breaking one of them.

- **Targeted changes.** Add, don't reformat. Keep diffs minimal and reviewable; don't restructure
  code you aren't changing.
- **Rules pass params explicitly.** A `.smk` rule reads from `config` and passes values through
  `params:`; the script reads `snakemake.params`, **never `config` directly**. After any change,
  cross-check that every param name matches between the `.smk` and its script.
- **Resources scale with retries.** Don't   hardcode resources in a rule. Every compute rule draws `mem_mb`/`runtime`/`threads` from
  `config["resources"]` (with a `default`), and `mem_mb` grows with the attempt number. 
- **The contract convention.** `obs["cell_id"]` must equal `obs_names` (as strings); downstream
  joins key on it. Head-produced contracts live at `SAMPLES_DIR/{sample}/{sample}_unfiltered.h5ad`
  (nested); decoupled contracts live at `contract_dir/{sample}.h5ad` (flat).
- **String indices when reading TSVs.** When a metadata/annotation TSV is read with
  `pd.read_csv(..., index_col=0)`, cast `df.index = df.index.astype(str)` before joining on
  `obs_names`. All-numeric cell ids are otherwise inferred as `int64` and silently misalign.
- **Guard plotting.** Wrap plot generation in `try/except` (Python) / `tryCatch` (R) so one failed
  figure doesn't crash the rule.
- **Palettes.** SpaceBlocks allows color palette customization directly from the config. Levels not referenced in config fall back to grey.
  When no palette is configured for a column, leave colours to scanpy and drop any stale `*_colors` from `uns`. 

## Adding things

**A new rule.** Add `workflow/rules/<name>.smk` (with `params:`) + `workflow/scripts/<name>.py`, a
`resources.<name>` block, and (if it introduces config keys) the corresponding entries in
`workflow/schemas/config.schema.yaml`. Wire it into the relevant named target
(`run_preprocessing` / `run_postprocessing` / `run_exploration`) in the `Snakefile`.

**A new platform (HeadBlock).** Add a `rules/*_<code>.smk` chain that emits the contract, a `<code>`
entry in `KNOWN_HEADS` and the `mode` schema enum, and its env under `envs/`. The CoreBlock needs no
changes. See the `xenium5k` head as a reference implementation.

## Validating a change

Run the applicable checks locally before pushing. The first two mirror CI and use the committed
`mode: decoupled` fixture under `.test/`:

```bash
# 1. Workflow parses + lints (against the committed test fixture)
snakemake -s workflow/Snakefile -d .test --lint --workflow-profile none

# 2. The decoupled CoreBlock DAG builds
snakemake -s workflow/Snakefile -d .test -n --workflow-profile none

# 3. Config still validates against the schema
python -c "import yaml,jsonschema; jsonschema.validate(yaml.safe_load(open('config/config.yaml')), yaml.safe_load(open('workflow/schemas/config.schema.yaml'))); print('Configuration schema validation passed')"

# 4. Docs build cleanly (only if you touched docs/)
pip install mkdocs-material pymdown-extensions
mkdocs build --strict
```

Before a catalogue release, generate the decoupled rule graph once manually:

```bash
snakemake -s workflow/Snakefile -d .test --forceall --rulegraph --workflow-profile none > pipeline_rulegraph.dot
```

Rendering the DOT file to SVG is optional and requires Graphviz:

```bash
dot -Tsvg pipeline_rulegraph.dot > images/rulegraph.svg
```

If you modify a HeadBlock, also dry-run the affected mode using a complete platform-specific configuration and representative mock or real inputs. Changing only `mode` in the generic `config/config.yaml` is not sufficient because its input paths are placeholders. Supported modes are `visiumhd`, `xenium5k`, `atera`, `merscope`, and `decoupled`.

Quick per-file sanity checks are cheap and worth it: `python -c "import ast; ast.parse(open('file.py').read())"`
for Python, and `yaml.safe_load` for any YAML you edit.

## Documentation

Docs live in `docs/` and are built with MkDocs Material; `mkdocs build --strict` **must pass**
(broken links and anchors fail the build). Inside `docs/`, use **relative** `.md` links: MkDocs
rewrites them and `--strict` validates them. In files rendered outside `docs/` (the root
`README.md`, `config/README.md`), use **absolute** `https://cbib.github.io/SpaceBlocks/...` URLs with
a trailing slash. `config/README.md` is also scraped by the Snakemake Workflow Catalog, so keep it
self-contained.

## Continuous integration

Two workflows run on every PR and must pass:

- **`tests`** — `snakemake --lint` and a decoupled dry-run against `.test/`.
- **`docs`** — `mkdocs build --strict`. Deployment to GitHub Pages is gated behind the repo variable
  `PUBLISH_DOCS=true` and only runs on `main`.

## Submitting a pull request

1. Branch from `main`, keep the change focused, and make sure the validation commands above pass.
2. Write a clear PR description: what changed and why. If behaviour changed, update the affected
   page under `docs/`.
3. For bug fixes, a one-line note of the root cause in the PR helps reviewers.

Questions or a larger design idea? Open an issue first so we can scope it together — especially
anything touching the contract or the HeadBlock/CoreBlock boundary. Thanks for contributing!
