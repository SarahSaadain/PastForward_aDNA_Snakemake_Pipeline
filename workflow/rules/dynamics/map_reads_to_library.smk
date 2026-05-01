
####################################################
# Snakemake rules
####################################################
_dyn_settings     = config.get("pipeline", {}).get("dynamics", {}).get("mapping", {}).get("settings", {})
_dyn_mapper       = _dyn_settings.get("mapper", "bwa-aln")
_BWA_ALN_DEFAULTS = "-n 0.01 -k 2 -l 1024 -o 2"  # Oliva et al. 2021 (10.1093/bib/bbab076)
_dyn_mapper_extra = _dyn_settings.get("mapper_extra_params", _BWA_ALN_DEFAULTS if _dyn_mapper == "bwa-aln" else "")

if _dyn_mapper == "minimap2":
    rule index_library_for_mapping_minimap2:
        input:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
        output:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.mmi",
        log:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg_minimap2_index.log"
        message: "Indexing SCG and Feature library {input} with minimap2"
        wrapper:
            "v9.3.0/bio/minimap2/index"

    rule map_reads_to_scg_feature_library_minimap2:
        input:
            query=["{species}/processed/reads/reads_merged/{individual}.fastq.gz"],
            target="{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.mmi",
        output:
            temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.bam"),
        log:
            "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_minimap2.log",
        message: "Mapping reads of {wildcards.individual} to {wildcards.species} SCG and Feature library with minimap2"
        params:
            extra="-ax sr " + _dyn_mapper_extra,
            sorting="none",
        threads: 10
        wrapper:
            "v9.3.0/bio/minimap2/aligner"

elif _dyn_mapper == "bwa-aln":
    rule index_library_for_mapping_bwa_aln:
        input:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
        output:
            multiext("{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        log:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg_bwa_aln_index.log"
        message: "Indexing SCG and Feature library {input} with BWA (for BWA ALN)"
        wrapper:
            "v9.3.0/bio/bwa/index"

    rule align_reads_to_library_bwa_aln:
        input:
            fastq="{species}/processed/reads/reads_merged/{individual}.fastq.gz",
            idx=multiext("{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        output:
            temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sai"),
        log:
            "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_bwa_aln.log",
        params:
            extra=_dyn_mapper_extra,
        threads: 10
        wrapper:
            "v9.3.0/bio/bwa/aln"

    rule map_reads_to_scg_feature_library_bwa_aln:
        input:
            fastq="{species}/processed/reads/reads_merged/{individual}.fastq.gz",
            sai="{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sai",
            idx=multiext("{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        output:
            temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.bam"),
        log:
            "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_bwa_samse.log"
        threads: 1
        wrapper:
            "v9.3.0/bio/bwa/samse"

else:
    # bwa-mem2 (default)
    rule index_library_for_mapping_bwa_mem2:
        input:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
        output:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.0123",
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.amb",
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.ann",
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.bwt.2bit.64",
            "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.pac",
        log:
            "{species}/processed/dynamics/lib/{feature_library}_and_scg_bwa_index.log"
        message: "Indexing SCG and Feature library {input} with BWA-MEM2"
        wrapper:
            "v9.3.0/bio/bwa-mem2/index"

    rule map_reads_to_scg_feature_library_bwa_mem2:
        input:
            reads=["{species}/processed/reads/reads_merged/{individual}.fastq.gz"],
            idx=multiext("{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        output:
            temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.bam"),
        log:
            "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_bwa.log",
        message: "Mapping reads of {wildcards.individual} to {wildcards.species} SCG and Feature library with BWA-MEM2"
        params:
            extra=_dyn_mapper_extra,
            sort="none",
            sort_order="coordinate",
            sort_extra="",
        threads: 10
        wrapper:
            "v9.3.0/bio/bwa-mem2/mem"


# Rule: Remove unmapped reads
rule remove_unmapped:
    # 2 Convert SAM to BAM
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.bam"
    output:
        bam=temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.no_unmapped.bam") # needs to be called .unsorted.bam otherwise snakemake had problems with unambigous names
    message: "Removing unmapped reads from BAM file for {input}"
    params:
        extra="-b -F 4",  # optional params string
    threads: 2
    wrapper:
        "v9.3.0/bio/samtools/view"

# Rule: Sort BAM file
rule  sort_bam_reads_to_library:
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.no_unmapped.bam"
    output:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam"
    message: "Sorting BAM file for {input}"
    log:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_sort_bam.log",
    threads: 8
    wrapper:
        "v9.3.0/bio/samtools/sort"

# Rule: Index BAM file
# SAMTOOLS doesn't parallelize the indexing work — it only parallelizes compression/decompression.
rule  index_bam_reads_to_library:
    # 4 Index BAM
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam"
    output:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam.bai"
    message: "Indexing BAM file for {input}"
    params:
        extra="",  # optional params string
    threads: 5
    wrapper:
        "v9.3.0/bio/samtools/index"
