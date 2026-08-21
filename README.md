# Snakemake version of the Galaxy workflow

This is a Snakemake translation of the Galaxy workflow stored in the accompanying `.ga` file. It processes paired-end sequencing data for reference-based viral variant calling and consensus generation, with support for multi-reference FASTA files and optional amplicon primer trimming.

## Pipeline overview

The pipeline runs the following steps in order:

| Step | Rule | Tool |
|---|---|---|
| Quality control & trimming | `qc_trim` | fastp |
| Reference indexing | `prepare_reference` | samtools, bwa |
| Read mapping | `map_reads` | bwa mem, samtools |
| Per-reference BAM splitting | `split_bam_by_reference` | samtools |
| Amplicon primer trimming (optional) | `trim_bam` | ivar trim |
| Consensus generation | `create_consensus` | samtools mpileup, ivar consensus |
| Variant calling | `call_variants_ivar` | samtools mpileup, ivar variants |
| ivar TSV → VCF conversion | `transform_ivar_to_vcf` | ivar_variants_to_vcf.py (Singularity) |
| Variant filtering | `filter_variants_lofreq` | lofreq filter |
| VCF normalization | `bcftools_normalize_variants` | bgzip, bcftools norm |

### Amplicon mode

If a `bed` file is provided for a sample, `ivar trim` is run after BAM splitting to remove amplicon primers. The `--ignore_strand_bias` flag is also applied automatically during VCF conversion.

### Multi-reference support

`reference_fasta` may contain multiple sequences. The pipeline:
1. Maps all reads against the full multi-reference FASTA
2. Splits the BAM by reference name
3. Runs all downstream steps (consensus, variant calling) independently per reference

## Outputs

All outputs are written under `res_dir` (configured in `config.yaml`), organized by sample:

```
<res_dir>/
  <sample>/
    qc/           # fastp output (trimmed reads, HTML/JSON report)
    mapping/      # BAMs per sample and per reference
    consensus/    # Per-reference consensus FASTA
    variants/     # Per-reference ivar TSV, VCF, filtered VCF, normalized VCF
    flags/        # Checkpoint flags
```

Logs are written to `logs/<sample>/` relative to the working directory.

## Sample sheet

The workflow is driven by a CSV sample sheet (path set in `config.yaml`).

**Required columns:**

| Column | Description |
|---|---|
| `sample` | Unique sample identifier |
| `reference_fasta` | Path to reference FASTA (single or multi-reference) |
| `r1` | Path to R1 FASTQ (gzipped) |
| `r2` | Path to R2 FASTQ (gzipped) |

**Optional columns:**

| Column | Description |
|---|---|
| `bed` | Path to BED file with amplicon primer coordinates |

Example:

| sample | reference_fasta | bed | r1 | r2 |
|---|---|---|---|---|
| sample1_trimmed | /data/refs/ref.fasta | /data/beds/primers.bed | /data/fastqs/sample1_R1.fastq.gz | /data/fastqs/sample1_R2.fastq.gz |
| sample1_not_trimmed | /data/refs/ref.fasta | | /data/fastqs/sample1_R1.fastq.gz | /data/fastqs/sample1_R2.fastq.gz |

## Configuration

All parameters are set in `config.yaml`:

| Parameter | Default | Description |
|---|---|---|
| `sample_sheet` | — | Path to sample sheet CSV |
| `res_dir` | — | Root output directory |
| `min_base_quality` | 20 | fastp minimum base quality (Phred) |
| `min_read_length` | 30 | fastp minimum read length |
| `fastp_extra` | `""` | Additional fastp arguments |
| `min_mapping_quality` | 20 | samtools view MAPQ filter |
| `bwa_mem_extra` | `""` | Additional bwa mem arguments |
| `min_length_after_trimming` | 50 | ivar trim minimum read length after primer removal |
| `min_variant_quality` | 30 | ivar minimum base quality for variant calling (Phred) |
| `min_allele_frequency` | 0.8 | ivar minimum allele frequency for variant calling |
| `min_allele_frequency_variant` | 0.1 | lofreq filter minimum allele frequency |
| `min_depth` | 20 | Minimum depth for variant calling and consensus |

## Conda environments

Each rule uses a dedicated conda environment defined under `requirements/`:

| File | Used by |
|---|---|
| `requirements_aln.yaml` | bwa, samtools, ivar |
| `requirements_mask_refs.yaml` | Biopython for `mask_refs.py` |
| `requirements_fastp.yaml` | fastp |
| `requirements_lofreq.yaml` | lofreq |
| `requirements_bcftools.yaml` | bcftools, bgzip |

The `transform_ivar_to_vcf` rule uses a Singularity container instead of conda (`biopython_matplotlib_pandas_python_pruned`).

## How to run

```bash
# Dry run
snakemake -n -p

# Standard run on CVG workstation:
# Replace <data_dir> with the path to your data directory (must contain fastqs, reference FASTAs, and BED files)
snakemake --cores 60 --use-conda --verbose --conda-frontend conda  --use-singularity --singularity-args "--bind <data_dir>:<data_dir>"

```

> **Note:** The explicit `--singularity-args --bind` is required on the CVG workstation because implicit directory mounting is not available. Set `<data_dir>` to the parent directory containing all your input data (fastqs, reference FASTAs, BED files). This path must match the directory paths used in `config.yaml` and the sample sheet.

## Workflow DAG

The DAG below shows an exemplary run using the sample sheet provided. The same sequencing data is processed twice — once with amplicon primer trimming and once without — by listing the sample under two different names and omitting the `bed` column for the untrimmed case. Both samples are mapped against a multi-FASTA reference containing the HSV-2 genes UL23 and UL30, illustrating the per-reference job expansion and the conditional primer trimming branch.

| sample | reference_fasta | bed | r1 | r2 |
|---|---|---|---|---|
| 1545554-HSV2_S82_trimmed | HSV_2_both.fasta | HSV_2_both.bed | 1545554-HSV2_S82_L001_R1_001.fastq.gz | 1545554-HSV2_S82_L001_R2_001.fastq.gz |
| 1545554-HSV2_S82_not_trimmed | HSV_2_both.fasta | | 1545554-HSV2_S82_L001_R1_001.fastq.gz | 1545554-HSV2_S82_L001_R2_001.fastq.gz |

![Workflow DAG](dag.svg)
