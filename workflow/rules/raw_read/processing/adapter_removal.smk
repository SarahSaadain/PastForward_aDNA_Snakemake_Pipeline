####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def remove_adapters_type_with_fastp_input_sample(wc):
    """
    Returns a list of read files for R1/R2 if available.
    If only R1 exists, returns a single-element list [R1].
    """

    return get_raw_reads_for_sample(wc.species, wc.sample)

def determine_reads_trimmed_final_input(wildcards):
    # Determine if the sample is paired-end or single-end
    reads = remove_adapters_type_with_fastp_input_sample(wildcards)
    if len(reads) == 2:
        # Paired-end: use the merged reads from fastp_pe
        return f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.pe.merged.fastq.gz"
    else:
        # Single-end: use the trimmed reads from fastp_se
        return f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.se.fastq.gz"
 
 

####################################################
# Snakemake rules
####################################################

# Rule: Adapter removal for single-end reads using fastp
rule remove_adapters_single_with_fastp:
    input:
        sample=remove_adapters_type_with_fastp_input_sample,
    output:
        trimmed=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.se.fastq.gz"),
        failed=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.se.failed.fastq.gz"),
        html="{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.se.html",
        json="{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.se.json",
    message: "Trimming adapters from single-end reads in {input.sample}"
    log:
        "{species}/logs/reads/reads_trimmed/{sample}_trimmed.se.log",
    params:
        adapters=lambda wc: (
            f"--adapter_sequence {config['pipeline']['raw_reads_processing']['adapter_removal']['settings']['adapters_sequences']['r1']}"
            if config['pipeline']['raw_reads_processing']['adapter_removal'].get('settings', {}).get('adapters_sequences', {}).get('r1')
            else ""
        ),  
        extra=lambda wc: (
            f"--length_required {config['pipeline']['raw_reads_processing']['adapter_removal'].get('settings', {}).get('min_length',0)} "
            f"--trim_poly_x 5 "
            f"--qualified_quality_phred {config['pipeline']['raw_reads_processing']['adapter_removal'].get('settings', {}).get('min_quality',0)} "
            f"--unqualified_percent_limit 40 "
            f"--n_base_limit 5"
        ),
    threads: 10
    wrapper:
        "v7.5.0/bio/fastp"
 
 
# Rule: Adapter removal for paired-end reads using fastp
rule remove_adapters_paired_with_fastp:
    input:
        sample=remove_adapters_type_with_fastp_input_sample,
    output:
        trimmed=[
            temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.R1.fastq.gz"),
            temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.R2.fastq.gz"),
        ],
        # Unpaired reads separately
        unpaired1=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.unpaired.R1.fastq.gz"),
        unpaired2=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.unpaired.R2.fastq.gz"),
        merged=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.fastq.gz"),
        failed=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.failed.fastq.gz"),
        html="{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.pe.html",
        json="{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.pe.json",
    message: "Trimming adapters from paired-end reads and merging for {input.sample}"
    log:
        "{species}/logs/reads/reads_trimmed/{sample}_trimmed.pe.log",
    params:
        adapters=lambda wc: (
            f"--adapter_sequence {config['pipeline']['raw_reads_processing']['adapter_removal']['settings']['adapters_sequences'].get('r1','')} "
            f"--adapter_sequence_r2 {config['pipeline']['raw_reads_processing']['adapter_removal']['settings']['adapters_sequences'].get('r2','')}"
            if config['pipeline']['raw_reads_processing']['adapter_removal'].get('settings', {}).get('adapters_sequences')
            else "--detect_adapter_for_pe"
        ),
        extra=lambda wc: (
            f"--length_required {config['pipeline']['raw_reads_processing']['adapter_removal'].get('settings', {}).get('min_length',0)} "
            f"--trim_poly_x 5 "
            f"--qualified_quality_phred {config['pipeline']['raw_reads_processing']['adapter_removal'].get('settings', {}).get('min_quality',0)} "
            f"--unqualified_percent_limit 40 "
            f"--n_base_limit 5 "
            f"--merge"
        ),   
    threads: 10
    wrapper:
        "v7.5.0/bio/fastp"

rule determine_reads_trimmed_final:
    input:
        determine_reads_trimmed_final_input,
    output:
        temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed_final.fastq.gz"),
    message: "Getting trimmed reads in {wildcards.sample}"
    shell:
        "mv {input} {output}"

rule merge_reads_trimmed_pe:
    input:
        merged="{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.fastq.gz",
        trimmed1="{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.R1.fastq.gz",
        trimmed2="{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.R2.fastq.gz",
        unpaired1="{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.unpaired.R1.fastq.gz",
        unpaired2="{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.unpaired.R2.fastq.gz",
    output:
        temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.merged.fastq.gz"),
    message: "Merging trimmed reads for paired-end"
    shell:
        "cat {input.merged} {input.trimmed1} {input.trimmed2} {input.unpaired1} {input.unpaired2} > {output}"