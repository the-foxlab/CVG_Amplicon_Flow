configfile: "config.yaml"

from pathlib import Path

import pandas as pd
from Bio import SeqIO


def _is_missing(value):
    text = str(value).strip()
    return not text or text.lower() == "nan"


def _read_fasta_names(fasta_path):
    names = []
    seen = set()
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            name = record.id
            if name in seen:
                raise ValueError(f"Duplicate FASTA record name '{name}' in {fasta_path}")
            seen.add(name)
            names.append(name)
    if not names:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    return names


def _read_bed_names(bed_path):
    names = []
    seen = set()
    with open(bed_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            if not name:
                continue
            if name not in seen:
                seen.add(name)
                names.append(name)
    return names


def _load_samples_from_sheet(sample_sheet):
    samples = {}
    table = pd.read_csv(sample_sheet, sep=',', engine="python")

    required_columns = {"sample", "reference_fasta", "r1", "r2"}
    if not required_columns.issubset(set(table.columns)):
        raise ValueError(
            "sample_sheet must contain columns: sample, reference_fasta, r1, r2"
        )

    has_bed_column = "bed" in set(table.columns)

    for row in table.itertuples(index=False):
        sample = str(row.sample).strip()
        if _is_missing(sample):
            continue

        reference_fasta = str(row.reference_fasta).strip()
        r1 = str(row.r1).strip()
        r2 = str(row.r2).strip()
        bed = str(row.bed).strip() if has_bed_column else ""

        if _is_missing(reference_fasta) or _is_missing(r1) or _is_missing(r2):
            raise ValueError(
                f"Sample '{sample}' has missing reference_fasta, r1 or r2 in sample_sheet"
            )

        if sample in samples:
            raise ValueError(f"Duplicate sample in sample_sheet: '{sample}'")

        reference_names = _read_fasta_names(reference_fasta)
        bed_names = []
        if not _is_missing(bed):
            bed_names = _read_bed_names(bed)
            unknown_bed_names = sorted(set(bed_names) - set(reference_names))
            if unknown_bed_names:
                raise ValueError(
                    f"BED file for sample '{sample}' contains references not present in FASTA: {unknown_bed_names}"
                )

        effective_references = bed_names if bed_names else reference_names
        if not effective_references:
            raise ValueError(
                f"Sample '{sample}' has no usable references after FASTA/BED parsing"
            )

        samples[sample] = {
            "reference_fasta": reference_fasta,
            "bed": "" if _is_missing(bed) else bed,
            "r1": r1,
            "r2": r2,
            "references": effective_references,
        }

    if not samples:
        raise ValueError("sample_sheet was provided but no samples were parsed")
    return samples


def _load_samples(cfg):
    sample_sheet = cfg.get("sample_sheet")
    if not sample_sheet:
        raise ValueError("sample_sheet is required")
    return _load_samples_from_sheet(sample_sheet)


def _reference_fasta(wildcards):
    return SAMPLES_MAP[wildcards.sample]["reference_fasta"]


def _masked_reference_fasta(wildcards):
    if SAMPLES_MAP[wildcards.sample]["bed"]:
        return _append_res_dir(f"{wildcards.sample}/references/{wildcards.sample}.masked.fasta")
    return SAMPLES_MAP[wildcards.sample]["reference_fasta"]


def _bed_file(wildcards):
    return SAMPLES_MAP[wildcards.sample].get("bed") or ""


def _append_res_dir(p):
    return str(Path(RES_DIR) / p)


RES_DIR = config.get('res_dir')
SAMPLES_MAP = _load_samples(config)
SAMPLES = sorted(SAMPLES_MAP.keys())
RUN_KEYS = [(sample, reference) for sample in SAMPLES for reference in SAMPLES_MAP[sample]["references"]]


rule all:
    input:
        [_append_res_dir(f"{sample}/consensus/{sample}_{reference}_consensus.fa") for sample, reference in RUN_KEYS],
        [_append_res_dir(f"{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.normalized.vcf") for sample, reference in RUN_KEYS[:]]


rule mask_reference:
    input:
        ref=_reference_fasta,
        bed=_bed_file
    output:
        ref=_append_res_dir("{sample}/references/{sample}.masked.fasta")
    log:
        stdout="logs/{sample}/mask_reference.stdout.log",
        stderr="logs/{sample}/mask_reference.stderr.log"
    params:
        script=str(Path(workflow.basedir) / "mask_refs.py")
    conda:
        "requirements/requirements_mask_refs.yaml"
    shell:
        """
        python {params.script} --fasta {input.ref} --bed {input.bed} --output {output.ref} \
            > {log.stdout} 2> {log.stderr}
        """


rule prepare_reference:
    input:
        ref=_masked_reference_fasta
    output:
        flag=_append_res_dir("{sample}/flags/reference_prepared.done")
    log:
        stdout="logs/{sample}/prepare_reference.stdout.log",
        stderr="logs/{sample}/prepare_reference.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools faidx {input.ref} > {log.stdout} 2> {log.stderr}
        bwa index {input.ref} > {log.stdout} 2>> {log.stderr}
        touch {output.flag}
        """

rule qc_trim:
    input:
        r1=lambda wildcards: SAMPLES_MAP[wildcards.sample]["r1"],
        r2=lambda wildcards: SAMPLES_MAP[wildcards.sample]["r2"]
    output:
        r1=_append_res_dir("{sample}/qc/{sample}_R1.fastq.gz"),
        r2=_append_res_dir("{sample}/qc/{sample}_R2.fastq.gz"),
        html=_append_res_dir("{sample}/qc/{sample}.html"),
        json=_append_res_dir("{sample}/qc/{sample}.json")
    log:
        stdout="logs/{sample}/qc_trim.stdout.log",
        stderr="logs/{sample}/qc_trim.stderr.log"
    params:
        min_quality=config.get("min_base_quality", 20),
        min_length=config.get("min_read_length", 30),
        extra=config.get("fastp_extra", "")
    threads:
        4
    conda:
        "requirements/requirements_fastp.yaml"
    shell:
        """
        fastp --thread {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} \
            --qualified_quality_phred {params.min_quality} \
            --length_required {params.min_length} --dont_eval_duplication \
            {params.extra} > {log.stdout} 2> {log.stderr}
        """

rule map_reads:
    input:
        r1=_append_res_dir("{sample}/qc/{sample}_R1.fastq.gz"),
        r2=_append_res_dir("{sample}/qc/{sample}_R2.fastq.gz"),
        ref_flag=_append_res_dir("{sample}/flags/reference_prepared.done"),
        ref=_masked_reference_fasta,
    output:
        sorted_bam=_append_res_dir("{sample}/mapping/{sample}.sorted.bam"),
        filtered_bam=_append_res_dir("{sample}/mapping/{sample}.sorted.filtered.bam")
    log:
        stdout="logs/{sample}/map_reads.stdout.log",
        stderr="logs/{sample}/map_reads.stderr.log"
    params:
        extra=config.get("bwa_mem_extra", ""),
        min_mapping_quality=config.get("min_mapping_quality", 30)
    threads:
        16
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        bwa mem -t {threads} -T 30 -h 5 {params.extra} {input.ref} {input.r1} {input.r2} 2> {log.stderr} \
            | samtools sort -o {output.sorted_bam} - > {log.stdout} 2>> {log.stderr}; 
        samtools view -b -f 3 -F 4 -q {params.min_mapping_quality} -o {output.filtered_bam} {output.sorted_bam} >> {log.stdout} 2>> {log.stderr};
        samtools index {output.filtered_bam} >> {log.stdout} 2>> {log.stderr};
        """


rule split_bam_by_reference:
    input:
        filtered_bam=_append_res_dir("{sample}/mapping/{sample}.sorted.filtered.bam")
    output:
        temp_bam=temp(_append_res_dir("{sample}/mapping/{sample}_{reference}.sorted.filtered.temp.bam")),
        bam=_append_res_dir("{sample}/mapping/{sample}_{reference}.sorted.filtered.bam"),
        bai=_append_res_dir("{sample}/mapping/{sample}_{reference}.sorted.filtered.bam.bai"),
    log: 
        stdout="logs/{sample}/split_bam_{reference}.stdout.log",
        stderr="logs/{sample}/split_bam_{reference}.stderr.log"
    params:
        exclude_refs = lambda wildcards: "|".join(
            ref for ref in SAMPLES_MAP[wildcards.sample]["references"] if ref != wildcards.reference
        ),
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools view -h {input.filtered_bam} {wildcards.reference} | egrep -v '{params.exclude_refs}' \
        | samtools view -bhS - > {output.temp_bam} 2> {log.stderr};
        samtools sort -o {output.bam} {output.temp_bam} >> {log.stdout} 2>> {log.stderr};
        samtools index {output.bam} >> {log.stdout} 2>> {log.stderr};
        """


rule trim_bam:
    input:
        split_bam=_append_res_dir("{sample}/mapping/{sample}_{reference}.sorted.filtered.bam"),
    output:
        temp_bam=temp(_append_res_dir("{sample}/mapping/{sample}_{reference}.temp.bam")),
        bam_trimmed=_append_res_dir("{sample}/mapping/{sample}_{reference}.sorted.filtered.trimmed.bam"),
        bai_trimmed=_append_res_dir("{sample}/mapping/{sample}_{reference}.sorted.filtered.trimmed.bam.bai")
    log: 
        stdout="logs/{sample}/trim_bam{reference}.stdout.log",
        stderr="logs/{sample}/trim_bam{reference}.stderr.log"
    params:
        bed_fn = lambda wildcards: SAMPLES_MAP[wildcards.sample]["bed"],
        prefix = lambda wildcards: _append_res_dir(f"{wildcards.sample}/mapping/{wildcards.sample}_{wildcards.reference}.temp.bam"),
        min_len_after_trimming = config['min_length_after_trimming']
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """        
        echo "$(date '+%Y-%m-%d %H:%M:%S') Starting ivar trim" >> {log.stdout};
        ivar trim -i {input.split_bam} -b {params.bed_fn} -p {params.prefix} \
        -x 0 -e -m {params.min_len_after_trimming} -q 20 -s 4  >> {log.stdout} 2>> {log.stderr};
        echo "$(date '+%Y-%m-%d %H:%M:%S') Done with ivar trim" >> {log.stdout};
        samtools sort -o {output.bam_trimmed} {output.temp_bam} >> {log.stdout} 2>> {log.stderr};
        samtools index {output.bam_trimmed} >> {log.stdout} 2>> {log.stderr};
        """


def get_bam_input(wildcards):
    if SAMPLES_MAP[wildcards.sample]["bed"]:
        bam_path = _append_res_dir(f"{wildcards.sample}/mapping/{wildcards.sample}_{wildcards.reference}.sorted.filtered.trimmed.bam")
    else:
        bam_path = _append_res_dir(f"{wildcards.sample}/mapping/{wildcards.sample}_{wildcards.reference}.sorted.filtered.bam")
    return bam_path


rule create_consensus:
    input:
        bam=get_bam_input
    output:
        _append_res_dir("{sample}/consensus/{sample}_{reference}_consensus.fa")
    params:
        min_quality=config.get("min_variant_quality", 20),
        min_freq=config.get("min_consensus_frequency", 0.8),
        min_depth=config.get("min_depth", 20),
        prefix = lambda wildcards: _append_res_dir(f"{wildcards.sample}/consensus/{wildcards.sample}_{wildcards.reference}_consensus")
    log:
        stdout="logs/{sample}/create_consensus_{reference}.stdout.log",
        stderr="logs/{sample}/create_consensus_{reference}.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools mpileup -A -a -d 0 -Q 0 {input.bam} 2> {log.stderr} | \
        ivar consensus -p {params.prefix} -q {params.min_quality} -t {params.min_freq} -c 0.75 -m {params.min_depth} -n N \
        > {log.stdout} 2>> {log.stderr};
        """


rule call_variants_ivar:
    input:
        bam=get_bam_input,
        ref=_masked_reference_fasta
    output:
        temp_fasta=_append_res_dir("{sample}/variants/{sample}_{reference}_temp_ref_file.fasta"),
        tsv=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.tsv")
    params:
        min_quality=config.get("min_variant_quality", 20),
        min_freq=config.get("min_allele_frequency", 0.8),
        min_depth=config.get("min_depth", 20),
        ivar_prefix=lambda wildcards: _append_res_dir(f"{wildcards.sample}/variants/{wildcards.sample}_{wildcards.reference}.ivar"),
        ref_name= lambda wildcards: wildcards.reference
    log:
        stdout="logs/{sample}/call_variants_ivar_{reference}.stdout.log",
        stderr="logs/{sample}/call_variants_ivar_{reference}.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools faidx {input.ref} {params.ref_name} > {output.temp_fasta} 2> {log.stderr};
        samtools mpileup -A -d 0 -B -Q 0 --reference {output.temp_fasta} {input.bam} 2>> {log.stderr} | \
            ivar variants -p {params.ivar_prefix} -q {params.min_quality} -t {params.min_freq} -m {params.min_depth} -r {input.ref} \
            > {log.stdout} 2>> {log.stderr};                
        """

rule transform_ivar_to_vcf:
    input:
        ivar_tsv=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.tsv"),
        temp_fasta=_append_res_dir("{sample}/variants/{sample}_{reference}_temp_ref_file.fasta"),
    output:        
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.vcf")
    log:
        stdout="logs/{sample}/transform_ivar_to_vcf_{reference}.stdout.log",
        stderr="logs/{sample}/transform_ivar_to_vcf_{reference}.stderr.log"
    container:
        'docker://community.wave.seqera.io/library/biopython_matplotlib_pandas_python_pruned:46d87e2ad1f8a063'
    params:
        ref_name= lambda wildcards: wildcards.reference,
        # --ignore_strand_bias for amplicon runs
        ignore_strand_bias = lambda wildcards: "--ignore_strand_bias" if SAMPLES_MAP[wildcards.sample]["bed"] else ""
    shell:
        """
        python external_scripts/ivar_variants_to_vcf.py {input.ivar_tsv} {output.vcf} --fasta {input.temp_fasta} \
        {params.ignore_strand_bias} > {log.stdout} 2> {log.stderr};
        """


rule filter_variants_lofreq:
    input:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.vcf")
    output:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.vcf")
    params:
        min_cov=config.get("min_coverage", 10),
        min_freq=config.get("min_allele_frequency_variant", 0.1)
    conda:
        "requirements/requirements_lofreq.yaml"
    log:
        stdout="logs/{sample}/filter_variants_lofreq_{reference}.stdout.log",
        stderr="logs/{sample}/filter_variants_lofreq_{reference}.stderr.log"
    shell:
        """
        lofreq filter -i {input.vcf} --no-defaults --verbose -v {params.min_cov} -V 0 -a 0.1 -A 0.0 -o {output.vcf} > {log.stdout} 2> {log.stderr};
        """

rule bcftools_normalize_variants:
    input:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.vcf"),
        ref=_masked_reference_fasta
    output:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.normalized.vcf")
    params:
        args="--check-ref w --site-win 1000 --sort pos --output-type 'v'"
    conda:
        "requirements/requirements_bcftools.yaml"
    threads:
        4
    log:
        stdout="logs/{sample}/filter_variants_lofreq_{reference}.stdout.log",
        stderr="logs/{sample}/filter_variants_lofreq_{reference}.stderr.log"
    shell:
        """
        bgzip -c {input.vcf} > {input.vcf}.gz 2> {log.stderr};
        bcftools index {input.vcf}.gz > {log.stdout} 2>> {log.stderr};
        bcftools norm --fasta-ref {input.ref} {params.args} \
        --threads {threads} {input.vcf}.gz > {output.vcf} 2>> {log.stderr};
        """




# rule de_novo_assembly:
#     input:
#         r1="results/{sample}/qc/{sample}_R1.fastq.gz",
#         r2="results/{sample}/qc/{sample}_R2.fastq.gz"
#     output:
#         "results/{sample}/assembly/{sample}_contigs.fasta"
#     params:
#         outdir="results/{sample}/assembly/spades"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         spades.py -1 {input.r1} -2 {input.r2} -o {params.outdir} --isolate
#         cp {params.outdir}/contigs.fasta {output}
#         """

# rule basic_metagenomics:
#     input:
#         r1="results/{sample}/qc/{sample}_R1.fastq.gz",
#         r2="results/{sample}/qc/{sample}_R2.fastq.gz"
#     output:
#         "results/{sample}/metagenomics/{sample}.txt"
#     shell:
#         """
#         mkdir -p results/{wildcards.sample}/metagenomics
#         echo "Placeholder for metagenomic profiling; replace with kraken2 or centrifuge" > {output}
#         """

# rule build_snpeff_db:
#     input:
#         ref="results/{sample}/mapping/{sample}.sorted.bam",
#         gff=GFF
#     output:
#         "results/snpeff/database.done"
#     shell:
#         """
#         mkdir -p results/snpeff
#         touch {output}
#         """
