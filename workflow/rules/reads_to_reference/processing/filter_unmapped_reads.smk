####################################################
# Snakemake rules for handling unmapped reads.
#
# These rules are inserted BETWEEN the last processing step
# (damage rescaling / deduplication / sorting) and get_final_bam,
# so the final BAM always reflects the chosen action.
#
# Statistics rules (samtools stats, qualimap) use get_pre_filter_bam_path
# (defined in get_final_bam.smk) so that unmapped-read metrics are
# always captured regardless of the chosen action.
#
# Config option (reference_processing.filter_unmapped_reads):
#   execute: false          # set to true to enable
#   action: "remove"        # "remove" | "extract_fastq" | "extract_fasta"
#
# action = "remove"
#   Produces a mapped-reads-only BAM as a temp intermediate, which
#   get_final_bam then copies to {individual}_{reference}_final.bam.
#   The original pre-filter BAM is left on disk (it is the last
#   processing-step output and is not marked temp).
#
# action = "extract_fastq" / "extract_fasta"
#   Extracts unmapped reads into a compressed FASTQ or FASTA file.
#   get_final_bam copies the unmodified pre-filter BAM to _final.bam.
####################################################

# ---------------------------------------------------------------------------
# action: "remove"
# Produce a mapped-reads-only BAM (flag -F 4 excludes unmapped reads).
# Marked temp because get_final_bam copies its content to _final.bam.
# ---------------------------------------------------------------------------
rule filter_mapped_only_bam:
    input:
        bam = _pre_filter_bam,
        bai = _pre_filter_bam + ".bai"
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_mapped_only.bam")
    params:
        extra="-F 4"   # exclude unmapped reads
    threads: 4
    log:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_mapped_only.log"
    message:
        "Filtering unmapped reads from BAM for {wildcards.individual} mapped to {wildcards.reference}."
    wrapper:
        "v9.3.0/bio/samtools/view"

rule index_mapped_only_bam:
    input:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_mapped_only.bam"
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_mapped_only.bam.bai")
    params:
        extra=""
    threads: 4
    log:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_mapped_only.bam.bai.log"
    message:
        "Indexing mapped-only BAM for {wildcards.individual} mapped to {wildcards.reference}."
    wrapper:
        "v9.3.0/bio/samtools/index"

# ---------------------------------------------------------------------------
# action: "extract_fastq"
# Extract unmapped reads (flag -f 4) into a compressed FASTQ file.
# The intermediate subset BAM is temp; the FASTQ is the kept output.
# ---------------------------------------------------------------------------
rule extract_unmapped_reads_subset_for_fastq:
    input:
        bam = _pre_filter_bam
    output:
        temp("{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_subset_fastq.bam")
    params:
        extra="-f 4"   # select only unmapped reads
    threads: 4
    log:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_subset_fastq.log"
    message:
        "Extracting unmapped read subset for FASTQ for {wildcards.individual} mapped to {wildcards.reference}."
    wrapper:
        "v9.3.0/bio/samtools/view"

rule convert_unmapped_reads_to_fastq:
    input:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_subset_fastq.bam"
    output:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped.fastq.gz"
    params:
        extra=""
    threads: 4
    log:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_fastq.log"
    message:
        "Converting unmapped reads to FASTQ for {wildcards.individual} mapped to {wildcards.reference}."
    wrapper:
        "v9.3.0/bio/samtools/fastq"

# ---------------------------------------------------------------------------
# action: "extract_fasta"
# Extract unmapped reads (flag -f 4) into a compressed FASTA file.
# ---------------------------------------------------------------------------
rule extract_unmapped_reads_subset_for_fasta:
    input:
        bam = _pre_filter_bam
    output:
        temp("{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_subset_fasta.bam")
    params:
        extra="-f 4"   # select only unmapped reads
    threads: 4
    log:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_subset_fasta.log"
    message:
        "Extracting unmapped read subset for FASTA for {wildcards.individual} mapped to {wildcards.reference}."
    wrapper:
        "v9.3.0/bio/samtools/view"

rule convert_unmapped_reads_to_fasta:
    input:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_subset_fasta.bam"
    output:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped.fasta.gz"
    params:
        extra=""
    threads: 4
    log:
        "{species}/processed/{reference}/unmapped/{individual}_{reference}_unmapped_fasta.log"
    message:
        "Converting unmapped reads to FASTA for {wildcards.individual} mapped to {wildcards.reference}."
    wrapper:
        "v9.3.0/bio/samtools/fasta"
