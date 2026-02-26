####################################################
# Snakemake rules
####################################################

# Rule: Quality filtering of reads using fastp
rule filter_reads_by_quality:
    input:
        sample = ["{species}/processed/reads/reads_trimmed/{sample}_trimmed_final.fastq.gz"]
    output:
        trimmed=temp("{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz"),
        failed=temp("{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.failed.fastq.gz"),
        html="{species}/results/reads/reads_quality_filtered/fastp_report/{sample}_quality_filtered.html",
        json="{species}/results/reads/reads_quality_filtered/fastp_report/{sample}_quality_filtered.json",
    message: "Quality filtering reads in {input.sample}"
    log:
        "{species}/logs/reads/reads_quality_filtered/{sample}_quality_filtered.log",
    params:
        extra=f"--disable_adapter_trimming --qualified_quality_phred {config['pipeline']['raw_reads_processing']['quality_filtering'].get('settings', {}).get('min_quality','15')} --length_required {config['pipeline']['raw_reads_processing']['quality_filtering'].get('settings', {}).get('min_length','15')} --unqualified_percent_limit 40 --n_base_limit 5"
    threads: 10
    wrapper:
        "v7.5.0/bio/fastp"