# Snakemake version of the Galaxy workflow

This directory contains a Snakemake translation of the Galaxy workflow stored in the accompanying .ga file.

## What is included
- Quality control with fastp
- Reference preparation with samtools + bwa
- Read mapping with bwa mem
- Variant calling with lofreq and ivar
- Consensus generation with ivar
- Optional de novo assembly with SPAdes
- Optional metagenomics placeholder

## How sample inputs work
The workflow is driven by a sample sheet with one row per sample.

Required sample sheet columns:
- sample
- reference_fasta
- r1
- r2

Optional sample sheet column:
- bed

Notes:
- Each sample appears only once in the sample sheet.
- `reference_fasta` can contain one reference sequence or multiple reference sequences.
- If `bed` is provided, it can also contain one or multiple reference-specific entries.
- Mapping first creates one BAM per sample, then splits it by reference name to produce per-reference BAMs (`sample` + `reference` wildcards).

## How to run
1. Adjust config.yaml:
	- set sample_sheet path
	- optionally set GFF
2. Install the required tools: fastp, bwa, samtools, lofreq, ivar, and optionally SPAdes.
3. Run:

```bash
snakemake -j 4
```

The workflow is intentionally structured as a scaffold that mirrors the Galaxy stages and can be customized further.
