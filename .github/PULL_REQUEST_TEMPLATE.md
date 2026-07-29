## What & why

<!-- One or two sentences: what changed and why. Link any related issue. -->

## Validation

<!-- Exact commands in CONTRIBUTING.md. Tick what applies. -->

- [ ] `snakemake --lint` passes (against the `.tests/` config)
- [ ] DAG builds — decoupled smoke test + any affected head mode (`snakemake -n`)
- [ ] Config still validates against `workflow/schemas/config.schema.yaml`
- [ ] `mkdocs build --strict` passes (if `docs/` was touched)
- [ ] Param names match between `.smk` and script; docs updated if behaviour changed
