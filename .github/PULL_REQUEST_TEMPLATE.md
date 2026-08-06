## What & why

<!-- One or two sentences: what changed and why. Link any related issue. -->

## Validation

<!-- Tick the checks that apply. Use the exact commands documented in CONTRIBUTING.md. -->

- [ ] Workflow lint passes against `.test/`
- [ ] Decoupled `.test/` dry-run builds the CoreBlock DAG
- [ ] Affected HeadBlock mode was dry-run with a complete mode-specific config and representative inputs, if applicable
- [ ] Config validates against `workflow/schemas/config.schema.yaml`
- [ ] `mkdocs build --strict` passes, if documentation was changed
- [ ] Catalogue rule graph was generated before a release, if applicable
- [ ] Parameter names match between `.smk` rules and scripts; documentation was updated if behaviour changed
