configfile: "config.yaml"

from pathlib import Path
import sys

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))

from snakefile_helpers import append_res_dir, build_bwa_files, load_samples, masked_reference_fasta


def _reference_fasta(wildcards):
    return SAMPLES_MAP[wildcards.sample]["reference_fasta"]


def _all_bwa_files(wildcards):
    return build_bwa_files(SAMPLES_MAP[wildcards.sample], RES_DIR)

    
def _bed_file(wildcards):
    return SAMPLES_MAP[wildcards.sample].get("bed") or ""


def _append_res_dir(p):
    return append_res_dir(RES_DIR, p)


def _masked_reference_fasta(wildcards):
    return masked_reference_fasta(SAMPLES_MAP[wildcards.sample], RES_DIR)


RES_DIR = config.get('res_dir')
SAMPLES_MAP = load_samples(config)
SAMPLES = sorted(SAMPLES_MAP.keys())
RUN_KEYS = [(sample, reference) for sample in SAMPLES for reference in SAMPLES_MAP[sample]["references"]]

# Check the README section "Hashing of references and bed files" for an explanation of this mechanism
HASH_to_files = {SAMPLES_MAP[sample]["reference_hash"]: [SAMPLES_MAP[sample].get("reference_fasta"), SAMPLES_MAP[sample].get("bed")] for sample in SAMPLES}

wildcard_constraints:
    ref_type="masked|unmasked"


rule all:
    input:
        [_append_res_dir(f"{sample}/consensus/{sample}_{reference}_consensus.fa") for sample, reference in RUN_KEYS],
        [_append_res_dir(f"{sample}/variants/{sample}_{reference_hash}.ivar.lofreq_filtered.normalized.vcf") for sample, reference_hash in RUN_KEYS[:]],
        [_append_res_dir(f"{sample}/visualization/{sample}_{reference}_plot.html") for sample, reference in RUN_KEYS]


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
        min_quality=config.get("fastp_min_base_quality", 20),
        min_length=config.get("fastp_min_read_length", 30),
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
    """
    Map reads to reference using BWA MEM, sort and filter the resulting BAM file.
    The filtering of samtools view is different to the galaxy pipeline. 
    This was adapted after discussions with Jonas Fuchs about the intended behavior. 
    """
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
        extra=config.get("bwa_mem_extra", ""),
        T=config.get("bwa_T", 30),
        h=config.get("bwa_h", 5),
        sam_f=config.get("samtools_filtering_post_mapping_f", "3"),
        sam_F=config.get("samtools_filtering_post_mapping_F", "4"),
        min_mapping_quality=config.get("min_mapping_quality", 30)
    threads:
        16
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        bwa mem -t {threads} -T {params.T} -h {params.h} {params.extra} {input.ref} {input.r1} {input.r2} 2> {log.stderr} \
            | samtools sort -o {output.sorted_bam} - > {log.stdout} 2>> {log.stderr}; 
        samtools view -b -f {params.sam_f} -F {params.sam_F} -q {params.min_mapping_quality} -o {output.filtered_bam} {output.sorted_bam} >> {log.stdout} 2>> {log.stderr};
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
        min_len_after_trimming = config['ivar_min_length_after_trimming'],
        trim_offset = config['ivar_trim_offset'],
        trim_q = config['ivar_trim_q'],
        trim_s = config['ivar_trim_s']
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """        
        echo "$(date '+%Y-%m-%d %H:%M:%S') Starting ivar trim" >> {log.stdout};
        ivar trim -i {input.split_bam} -b {params.bed_fn} -p {params.prefix} \
        -x {params.trim_offset} -e -m {params.min_len_after_trimming} -q {params.trim_q} -s {params.trim_s}  >> {log.stdout} 2>> {log.stderr};
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
        min_quality=config.get("cons_ivar_min_variant_quality", 30),
        min_freq=config.get("cons_ivar_min_consensus_frequency", 0.75),
        min_insertion_freq=config.get("cons_ivar_min_insertion_frequency", 0.75),
        min_depth=config.get("cons_ivar_min_depth", 20),
        prefix = lambda wildcards: _append_res_dir(f"{wildcards.sample}/consensus/{wildcards.sample}_{wildcards.reference}_consensus"),
        mpileup_A = "-A" if config.get("cons_mpileup_A", True) else "",
        mpileup_a = "-a" if config.get("cons_mpileup_a", True) else "",
        mpileup_d = config.get("cons_mpileup_extra", 0),
        mpileup_Q = config.get("cons_mpileup_Q", 0)
    log:
        stdout="logs/{sample}/create_consensus_{reference}.stdout.log",
        stderr="logs/{sample}/create_consensus_{reference}.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools mpileup {params.mpileup_A} {params.mpileup_a} -d {params.mpileup_d} -Q {params.mpileup_Q} \
         {input.bam} 2> {log.stderr} | \
        ivar consensus -p {params.prefix} -q {params.min_quality} -t {params.min_freq} -c {params.min_insertion_freq} \
         -m {params.min_depth} -n N > {log.stdout} 2>> {log.stderr};
        """

rule call_variants_ivar:
    input:
        bam=get_bam_input,
        ref=_masked_reference_fasta
    output:
        temp_fasta=_append_res_dir("{sample}/variants/{sample}_{reference}_temp_ref_file.fasta"),
        tsv=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.tsv")
    params:
        mpileup_A = "-A" if config.get("ivar_variant_mpileup_A", True) else "",
        mpileup_d = config.get("ivar_variant_mpileup_d", 0),
        mpileup_B = "-B" if config.get("ivar_variant_mpileup_B", True) else "",
        mpileup_Q = config.get("ivar_variant_mpileup_Q", 0),
        ivar_prefix=lambda wildcards: _append_res_dir(f"{wildcards.sample}/variants/{wildcards.sample}_{wildcards.reference}.ivar"),

        ivar_q=config.get("ivar_variant_min_variant_quality", 30),
        ivar_t=config.get("ivar_variant_t", 0.0),
        ref_name= lambda wildcards: wildcards.reference
    log:
        stdout="logs/{sample}/call_variants_ivar_{reference}.stdout.log",
        stderr="logs/{sample}/call_variants_ivar_{reference}.stderr.log"
    conda:
        'requirements/requirements_aln.yaml'
    shell:
        """
        samtools faidx {input.ref} {params.ref_name} > {output.temp_fasta} 2> {log.stderr};
        samtools mpileup {params.mpileup_A} -d {params.mpileup_d} {params.mpileup_B} -Q {params.mpileup_Q} \
        --reference {output.temp_fasta} {input.bam} 2>> {log.stderr} | \
        ivar variants -p {params.ivar_prefix} -q {params.ivar_q} -t {params.ivar_t} -r {input.ref} \
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
        lofreq_no_defaults = "--no-defaults " if config.get("lofreq_filter_no_defaults", True) else "",
        lofreq_v = config.get("lofreq_filter_v", 20),
        lofreq_V = config.get("lofreq_filter_V", 0),
        lofreq_a = config.get("lofreq_filter_a", 0.1),
        lofreq_A = config.get("lofreq_filter_A", 0.0)
    conda:
        "requirements/requirements_lofreq.yaml"
    log:
        stdout="logs/{sample}/filter_variants_lofreq_{reference}.stdout.log",
        stderr="logs/{sample}/filter_variants_lofreq_{reference}.stderr.log"
    shell:
        """
        lofreq filter -i {input.vcf} {params.lofreq_no_defaults} --verbose -v {params.lofreq_v} -V {params.lofreq_V} \
            -a {params.lofreq_a} -A {params.lofreq_A} -o {output.vcf} > {log.stdout} 2> {log.stderr};
        """

rule bcftools_normalize_variants:
    input:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.vcf"),
        ref=_masked_reference_fasta
    output:
        vcf=_append_res_dir("{sample}/variants/{sample}_{reference}.ivar.lofreq_filtered.normalized.vcf")
    params:
        check_ref=config.get("bcftools_normalize_check_ref", "w"),
        site_win=config.get("bcftools_normalize_site_win", 1000),
        sort=config.get("bcftools_normalize_sort", "pos"),
        out_type=config.get("bcftools_normalize_out_type", "v")
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
        bcftools norm --fasta-ref {input.ref} --check-ref {params.check_ref} --site-win {params.site_win} \
        --sort {params.sort} --output-type {params.out_type} \
        --threads {threads} {input.vcf}.gz > {output.vcf} 2>> {log.stderr};
        """

rule visualize_bam_bamdash:
    input:
        bam=get_bam_input,
        bai=lambda wildcards: f"{get_bam_input(wildcards)}.bai",
    output:
        html=_append_res_dir("{sample}/visualization/{sample}_{reference}_plot.html"),
    log:
        stdout=str(Path(workflow.basedir) / "logs/{sample}/visualize_bam_bamdash_{reference}.stdout.log"),
        stderr=str(Path(workflow.basedir) / "logs/{sample}/visualize_bam_bamdash_{reference}.stderr.log")
    params:
        ref_id = lambda wildcards: wildcards.reference,
        prefix = _append_res_dir("{sample}/visualization/{sample}_{reference}_plot"),
        bs = config.get("bamdash_bs", 10)
    conda:
        "requirements/requirements_bamdash.yaml"
    shell:
        """
        bamdash --bam {input.bam} -r {params.ref_id} -bs {params.bs} -p {params.prefix} >> {log.stdout} 2>> {log.stderr};
        """