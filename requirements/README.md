# Conda requirements from Galaxy workflow metadata

Version provenance:
- Exact values extracted recursively from the Galaxy workflow are recorded in:
  - ../galaxy_tool_versions.tsv
  - ../galaxy_tool_versions_unique.tsv
- Those files contain exact Galaxy tool wrapper versions (`tool_version` / `tool_id`).

Conda files in this directory:
- requirements-core.yaml
- requirements-reference-annotation.yaml
- requirements-denovo.yaml
- requirements-metagenomics-kraken2-2.17.1.yaml
- requirements-metagenomics-kraken2-2.1.1.yaml

Why multiple files:
- The Galaxy workflow uses conflicting versions of some wrappers across different subworkflows (notably Kraken2 wrapper versions `2.17.1+galaxy0` and `2.1.1+galaxy1`).
- Separate envs keep those constraints isolated.

Important note:
- Galaxy wrapper versions are exact from the workflow export.
- Conda package pins are mapped to corresponding software names and may need adjustment if a specific build is unavailable on your platform.

Create an environment:

```bash
conda env create -f requirements/requirements-core.yaml
```
