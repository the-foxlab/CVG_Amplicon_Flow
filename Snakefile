configfile: "config.yaml"

from pathlib import Path
from galaxy_port.foundation import PipelineConfigResolver, build_pipeline_samples_map


def _load_samples(cfg):
    sample_sheet = cfg.get("sample_sheet")
    if not sample_sheet:
        raise ValueError("sample_sheet is required")
    return build_pipeline_samples_map(sample_sheet)


def _reference_fasta(wildcards):
    return SAMPLES_MAP[wildcards.sample]["reference_fasta"]


def _all_bwa_files(wildcards):
    sample = SAMPLES_MAP[wildcards.sample]
    hash_value = sample["reference_hash"]

    ref_type = "masked" if sample["bed"] else "unmasked"

    ref = _append_res_dir(
        f"references/{hash_value}/{hash_value}.{ref_type}.fasta"
    )

    return {
        "ref": ref,
        "ref_fai": f"{ref}.fai",
        "ref_amb": f"{ref}.amb",
        "ref_ann": f"{ref}.ann",
        "ref_bwt": f"{ref}.bwt",
        "ref_pac": f"{ref}.pac",
        "ref_sa": f"{ref}.sa",
    }

    
def _bed_file(wildcards):
    return SAMPLES_MAP[wildcards.sample].get("bed") or ""
def _append_res_dir(p):
    return str(Path(RES_DIR) / p)

def _masked_reference_fasta(wildcards):
    sample = SAMPLES_MAP[wildcards.sample]
    hash_value = sample["reference_hash"]
    ref_type = "masked" if sample["bed"] else "unmasked"
    return _append_res_dir(f"references/{hash_value}/{hash_value}.{ref_type}.fasta")


def _tool_args(wildcards, tool_id, legacy_extra_key=None):
    return CONFIG_RESOLVER.render_tool_args(wildcards.sample, tool_id, legacy_extra_key)


RES_DIR = config.get('res_dir')
CONFIG_RESOLVER = PipelineConfigResolver(config)
SAMPLES_MAP = _load_samples(config)
SAMPLES = sorted(SAMPLES_MAP.keys())
RUN_KEYS = [(sample, reference) for sample in SAMPLES for reference in SAMPLES_MAP[sample]["references"]]
HASH_to_files = {SAMPLES_MAP[sample]["reference_hash"]: [SAMPLES_MAP[sample].get("reference_fasta"), 
                                                        SAMPLES_MAP[sample].get("bed")] for sample in SAMPLES}

wildcard_constraints:
    ref_type="masked|unmasked"


rule all:
    input:
        [_append_res_dir(f"{sample}/consensus/{sample}_{reference}_consensus.fa") for sample, reference in RUN_KEYS],
        [_append_res_dir(f"{sample}/variants/{sample}_{reference_hash}.ivar.lofreq_filtered.normalized.vcf") for sample, reference_hash in RUN_KEYS[:]],
        [_append_res_dir(f"{sample}/visualization/{sample}_{reference}/{sample}_{reference}.html") for sample, reference in RUN_KEYS]

def _get_masking_input(wildcards):
    [fasta_fn, bed_fn] = HASH_to_files[wildcards.hash_value]
    assert bed_fn, (
        f"Reference {fasta_fn} has no BED file, "
        "but mask_reference was triggered"
    )
    return { "ref": fasta_fn, "bed": bed_fn }

rule mask_reference:
    input:
        unpack(_get_masking_input)
    output:
        ref = _append_res_dir("references/{hash_value}/{hash_value}.masked.fasta")
    log:
        stdout="logs/mask_reference/{hash_value}.stdout.log",
        stderr="logs/mask_reference/{hash_value}.stderr.log"
    params:
        script=str(Path(workflow.basedir) / "scripts/mask_refs.py")
    conda:
        "requirements/requirements_mask_refs.yaml"
    shell:
        """
        python {params.script} --fasta {input.ref} --bed {input.bed} --output {output.ref} \
            > {log.stdout} 2> {log.stderr}
        """

def _get_reference_fasta(wildcards):
    [fasta_fn, bed_fn] = HASH_to_files[wildcards.hash_value]
    assert not bed_fn, f"Reference {wildcards.hash_value} has a BED file, but the rule mask_reference is not triggered"
    return fasta_fn

rule copy_reference:
    input:
        ref=_get_reference_fasta,
    output:
        ref = _append_res_dir("references/{hash_value}/{hash_value}.unmasked.fasta")
    log:
        stdout="logs/copy_reference/{hash_value}.stdout.log",
        stderr="logs/copy_reference/{hash_value}.stderr.log"
    shell:
        """
        cp {input.ref} {output.ref} > {log.stdout} 2> {log.stderr}
        """


rule prepare_reference:
    input:
        ref = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta")
    output:
        ref_fai = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta.fai"),
        ref_amb = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta.amb"),
        ref_ann = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta.ann"),
        ref_bwt = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta.bwt"),
        ref_pac = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta.pac"),
        ref_sa = _append_res_dir("references/{hash_value}/{hash_value}.{ref_type}.fasta.sa"),
    log:
        stdout="logs/prepare_reference/{hash_value}.{ref_type}.stdout.log",
        stderr="logs/prepare_reference/{hash_value}.{ref_type}.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools faidx {input.ref} > {log.stdout} 2> {log.stderr}
        bwa index {input.ref} >> {log.stdout} 2>> {log.stderr}
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
        tool_args=lambda wildcards: _tool_args(wildcards, "fastp", "fastp_extra")
    threads:
        4
    conda:
        "requirements/requirements_fastp.yaml"
    shell:
        """
        fastp --thread {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} \
            {params.tool_args} > {log.stdout} 2> {log.stderr}
        """

rule map_reads:
    input:
        unpack(_all_bwa_files),
        r1=_append_res_dir("{sample}/qc/{sample}_R1.fastq.gz"),
        r2=_append_res_dir("{sample}/qc/{sample}_R2.fastq.gz"),
    output:
        sorted_bam=_append_res_dir("{sample}/mapping/{sample}.sorted.bam"),
        filtered_bam=_append_res_dir("{sample}/mapping/{sample}.sorted.filtered.bam"),
        filtered_bai=_append_res_dir("{sample}/mapping/{sample}.sorted.filtered.bam.bai")
    log:
        stdout="logs/{sample}/map_reads.stdout.log",
        stderr="logs/{sample}/map_reads.stderr.log"
    params:
        bwa_mem_args=lambda wildcards: _tool_args(wildcards, "bwa_mem", "bwa_mem_extra"),
        samtools_filter_args=lambda wildcards: _tool_args(wildcards, "samtools_view_filter")
    threads:
        16
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        bwa mem -t {threads} {params.bwa_mem_args} {input.ref} {input.r1} {input.r2} 2> {log.stderr} \
            | samtools sort -o {output.sorted_bam} - > {log.stdout} 2>> {log.stderr}; 
        samtools view -b {params.samtools_filter_args} -o {output.filtered_bam} {output.sorted_bam} >> {log.stdout} 2>> {log.stderr};
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
        tool_args = lambda wildcards: _tool_args(wildcards, "ivar_trim")
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """        
        echo "$(date '+%Y-%m-%d %H:%M:%S') Starting ivar trim" >> {log.stdout};
        ivar trim -i {input.split_bam} -b {params.bed_fn} -p {params.prefix} \
        {params.tool_args}  >> {log.stdout} 2>> {log.stderr};
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
        tool_args=lambda wildcards: _tool_args(wildcards, "ivar_consensus"),
        prefix = lambda wildcards: _append_res_dir(f"{wildcards.sample}/consensus/{wildcards.sample}_{wildcards.reference}_consensus")
    log:
        stdout="logs/{sample}/create_consensus_{reference}.stdout.log",
        stderr="logs/{sample}/create_consensus_{reference}.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools mpileup -A -a -d 0 -Q 0 {input.bam} 2> {log.stderr} | \
        ivar consensus -p {params.prefix} {params.tool_args} \
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
        tool_args=lambda wildcards: _tool_args(wildcards, "ivar_variants"),
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
            ivar variants -p {params.ivar_prefix} {params.tool_args} -r {input.ref} \
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
        # --ignore_strand_bias for amplicon runs, for now i always ignore it
        options = "--ignore_strand_bias"
    shell:
        """
        python external_scripts/ivar_variants_to_vcf.py {input.ivar_tsv} {output.vcf} --fasta {input.temp_fasta} \
        {params.options} > {log.stdout} 2> {log.stderr};
        """


rule filter_variants_lofreq:
    input:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.vcf")
    output:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.vcf")
    params:
        tool_args=lambda wildcards: _tool_args(wildcards, "lofreq_filter")
    conda:
        "requirements/requirements_lofreq.yaml"
    log:
        stdout="logs/{sample}/filter_variants_lofreq_{reference}.stdout.log",
        stderr="logs/{sample}/filter_variants_lofreq_{reference}.stderr.log"
    shell:
        """
        lofreq filter -i {input.vcf} {params.tool_args} -o {output.vcf} > {log.stdout} 2> {log.stderr};
        """

rule bcftools_normalize_variants:
    input:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.vcf"),
        ref=_masked_reference_fasta
    output:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.normalized.vcf")
    params:
        args=lambda wildcards: _tool_args(wildcards, "bcftools_norm")
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

rule visualize_bam_bamdash:
    input:
        bam=get_bam_input,
        bai=lambda wildcards: f"{get_bam_input(wildcards)}.bai",
    output:
        html=_append_res_dir("{sample}/visualization/{sample}_{reference}/{sample}_{reference}.html")
    log:
        stdout=str(Path(workflow.basedir) / "logs/{sample}/visualize_bam_bamdash_{reference}.stdout.log"),
        stderr=str(Path(workflow.basedir) / "logs/{sample}/visualize_bam_bamdash_{reference}.stderr.log")
    params:
        ref_id = lambda wildcards: wildcards.reference,
        out_dir = _append_res_dir("{sample}/visualization/{sample}_{reference}"),
        tool_args = lambda wildcards: _tool_args(wildcards, "bamdash")
    conda:
        "requirements/requirements_bamdash.yaml"
    shell:
        """
        mkdir -p {params.out_dir} > {log.stdout} 2> {log.stderr};
        cd {params.out_dir} >> {log.stdout} 2>> {log.stderr};
        bamdash --bam {input.bam} -r {params.ref_id} {params.tool_args}  >> {log.stdout} 2>> {log.stderr};
        mv {params.ref_id}_plot.html {output.html} >> {log.stdout} 2>> {log.stderr}
        """