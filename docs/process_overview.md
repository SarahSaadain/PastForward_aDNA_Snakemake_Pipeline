# Past Forward aDNA Pipeline — Process Overview

This document describes the processing logic of the Past Forward aDNA Pipeline across its three main modules: raw read processing, reference processing, and dynamics processing. A fourth module handles summary reporting across all results. All pipeline behaviour is controlled through `config/config.yaml`.

![Pipeline Overview](img/aDNA_pipeline_process.svg)

---

## Module 1 — Raw Read Processing

This module takes raw sequencing data and prepares it for all downstream analyses. It runs on a per-sample basis and produces clean, merged reads ready for mapping.

### Prepare Raw Reads

Before any processing begins, raw FASTQ files are moved into a standardised directory structure. This staging step ensures that all subsequent rules consistently locate their input data regardless of where files were originally placed on disk.

### Adapter Removal

Adapter sequences are removed from raw reads using **fastp**. This step can be turned on or off with **`pipeline > raw_reads_processing > adapter_removal > execute`**. The pipeline automatically detects whether a sample is single-end or paired-end by checking for the presence of an R2 file, and each mode is handled by a dedicated rule.

For single-end data, an adapter sequence can optionally be provided via **`pipeline > raw_reads_processing > adapter_removal > settings > adapters_sequences > r1`**. If left empty, fastp performs automatic adapter detection. For paired-end data, separate sequences can be set for R1 and R2 via **`adapters_sequences > r1`** and **`adapters_sequences > r2`**; if neither is provided, fastp's built-in paired-end detection is used instead. Beyond adapters, the step enforces a minimum read length (**`settings > min_length`**, default 0) and a minimum per-base quality score (**`settings > min_quality`**, default 0). Poly-X tail trimming is always active.

Paired-end reads receive additional treatment: overlapping read pairs are merged into single reads, which is particularly important for ancient DNA where fragment lengths are often shorter than the combined read length. The final output for PE samples concatenates the merged reads, the unmerged trimmed R1 and R2, and any unpaired reads into a single FASTQ file — ensuring no data is discarded. Fastp generates an HTML and JSON report for every sample, which feed into later QC aggregation.

### Quality Filtering

After adapter removal, reads go through a second dedicated quality filtering step, again using **fastp**, controlled by **`pipeline > raw_reads_processing > quality_filtering > execute`**. Adapter trimming is explicitly disabled here so the step acts purely as a quality gate. The minimum quality score (**`pipeline > raw_reads_processing > quality_filtering > settings > min_quality`**, default 15) and minimum read length (**`settings > min_length`**, default 15) are configured independently from the adapter removal thresholds, allowing for separate tuning of each stage. The step produces its own HTML and JSON reports per sample.

### Merge by Individual

Many sequencing projects split a single individual across multiple sequencing runs or lanes, each producing a separate FASTQ file. This step concatenates all quality-filtered samples belonging to the same individual into a single merged FASTQ. Individual identity is derived from the sample filename — everything before the first underscore is treated as the individual identifier. The merged file serves as the single input to all downstream reference mapping and dynamics analyses.

### Read Count Statistics

Read counts are tracked at each processing stage — raw, after adapter removal, and after quality filtering — for every sample. These per-sample counts are collected into a species-level summary CSV. The counts are later used for QC plotting and are incorporated into the MultiQC summary reports.

### Quality Checks with FastQC and MultiQC

**FastQC** is run at four points in the workflow to track how read quality evolves through processing: on the raw reads, on the adapter-trimmed reads, on the quality-filtered reads, and on the final merged per-individual reads.

After FastQC, **MultiQC** aggregates the results from each stage into a single HTML report per stage per species, giving a consolidated quality overview for all samples at a glance.

### Read Count Plots

Two R-based plots are generated from the read count statistics: one showing total read counts per sample across all processing stages, and one comparing read counts grouped by individual. These plots provide a quick visual check of how much data was retained through the processing steps.

### Contamination Analysis

The Past Forward aDNA Pipeline supports two contamination detection tools, both operating on the quality-filtered reads. The entire contamination analysis block is controlled by **`pipeline > raw_reads_processing > contamination_analysis > execute`**, and each tool can additionally be toggled individually.

**eCMSD** assesses contamination by mapping reads against a curated mitochondrial reference database and reporting the proportional contribution of different taxa. It is enabled via **`tools > ecmsd > execute`**. The path to the eCMSD executable is set with **`tools > ecmsd > settings > executable`** and the conda environment can be customised via **`settings > conda_env`**. Several analysis parameters are configurable: **`settings > Binsize`** (default 1000) controls the read binning resolution, **`settings > RMUS_threshold`** (default 0.15) filters out marginal hits, **`settings > mapping_quality`** (default 20) sets the minimum mapping quality for a read to be considered, and **`settings > taxonomic_hierarchy`** (default `genus`) determines at which taxonomic level results are reported. Results from all samples belonging to the same individual are merged into a combined summary file.

**Centrifuge** performs k-mer-based taxonomic classification against a user-provided database and is enabled via **`tools > centrifuge > execute`**. The **`settings > index`** parameter is required and must point to the full path of the Centrifuge database index prefix. The conda environment can be overridden with **`settings > conda_env`**. Beyond raw classification, the pipeline derives proportional taxon abundance and extracts the top 10 taxa ranked by both total and unique read assignments.

---

## Module 2 — Reference Processing

This module maps reads to one or more reference genomes and produces a final analysis-ready BAM for each individual, along with a comprehensive suite of mapping statistics and coverage metrics. All steps run per individual per reference. The entire module is gated by **`pipeline > reference_processing > execute`** (default `true`).

### Reference Preparation

Before mapping can begin, the reference genome is standardised to a `.fa` extension (regardless of whether the input is `.fna`, `.fasta`, or `.fa`) and two index files are created: a BWA index for read mapping and a samtools FAI index used by deduplication and damage analysis. Multiple reference genomes per species are supported — each is given a unique identifier derived from its filename and processed entirely independently. BWA indexing results are cached so the index is only rebuilt if the reference changes.

### Read Mapping

The merged per-individual reads are mapped to the reference using **BWA mem**. The resulting alignments are immediately sorted by coordinate and indexed. The unsorted BAM is discarded to save disk space. At this point the sorted BAM represents all mapped reads including duplicates, and is used as-is for library complexity estimation before any duplicate removal or damage rescaling.

### Deduplication

PCR and sequencing duplicates are removed using **DeDup**, a tool specifically designed for ancient DNA that correctly handles merged single-stranded reads. Deduplication is controlled by **`pipeline > reference_processing > deduplication > execute`** (default `true`). If disabled, the sorted BAM from mapping is passed directly to subsequent steps.

Because DeDup can be memory-intensive on reference genomes with many contigs, the pipeline uses a divide-and-conquer approach: contigs are grouped into clusters, the sorted BAM is split by cluster, each cluster is deduplicated independently, and the results are merged back together. The cluster size range is configurable via **`deduplication > settings > min_contigs_per_cluster`** (default 10) and **`settings > max_contigs_per_cluster`** (default 500). Each deduplication run produces a histogram and a JSON statistics file; the per-cluster JSON files are merged into a single summary that feeds into downstream QC reporting.

### DNA Damage Analysis and BAM Rescaling

Ancient DNA is characterised by cytosine deamination, appearing as C→T substitutions at the 5' end and G→A substitutions at the 3' end of reads. This step profiles those damage patterns and optionally corrects for them. It is controlled by **`pipeline > reference_processing > damage_rescaling > execute`** (default `true`).

The input BAM for this step is selected dynamically: if deduplication was enabled the deduplicated BAM is used, otherwise the sorted BAM from mapping. The tool **mapDamage2** is run on the selected BAM to estimate damage patterns and rescale base quality scores accordingly. The rescaled BAM is then sorted and indexed, and the unsorted rescaled BAM is discarded. The mapDamage2 output directory, including the rescaled BAM and all statistics files, is copied into the summary folder structure to be included in the MultiQC report.

The input BAM for this step is selected dynamically: if deduplication was enabled the deduplicated BAM is used, otherwise the sorted BAM from mapping.

### Final BAM

After all optional processing steps, a single canonical `_final.bam` is produced by selecting the most-processed available BAM following a priority chain: rescaled BAM (if **`damage_rescaling > execute`** is true) → deduplicated BAM (if **`deduplication > execute`** is true) → sorted BAM. This ensures all downstream analytics always work from a consistently named file regardless of which steps were enabled.

### Coverage Analysis

Coverage analysis is controlled by **`pipeline > reference_processing > coverage_analysis > execute`** (default `true`). The final BAM is interrogated using **samtools depth**, which reports per-position depth across the entire reference including zero-coverage sites. A custom Python script then analyses the depth output to compute coverage breadth at multiple depth thresholds and mean depth statistics per individual. All individual results are combined into species-level summary tables, and the data is reformatted into MultiQC-compatible custom content files.

### Mapping Statistics

Four additional tools are run to produce complementary QC metrics:

**samtools stats** generates a comprehensive statistics file from the final BAM, used among other things to extract the fraction of reads that successfully mapped to the reference (the endogenous content).

**Qualimap bamqc** produces an HTML quality report covering mapping rates, coverage uniformity, GC content, and insert size distributions.

**Picard MarkDuplicates** is run on the sorted (pre-dedup) BAM to provide duplicate metrics for QC comparison purposes.

**Preseq** estimates library complexity and extrapolates sequencing yield. Critically, preseq must be run on the sorted BAM *before* deduplication and damage rescaling — running it afterwards removes the duplication information that preseq requires, leading to invalid estimates. Both a complexity curve and an lc_extrap extrapolation are generated.

### Endogenous Read Fraction

The samtools stats output is parsed to determine how many reads mapped to the reference genome. This endogenous content is reported per individual and combined across all individuals into a species-level file. Together with the raw read counts from Module 1, this provides a full picture of how much of the sequenced data is target-derived.

### Plots

A set of R-based plots summarise the mapping results visually. Coverage breadth and depth are each shown as both a violin plot and a bar chart, comparing all individuals for a given species and reference. Endogenous read fractions are shown as a bar chart per individual, and a combined plot shows raw read counts alongside endogenous fractions to provide full context from raw sequencing through to usable data.

### MultiQC BAM Report

All QC outputs for a given individual and reference are aggregated into a single **MultiQC** HTML report. This includes the fastp reports from adapter removal and quality filtering, FastQC of merged reads, contamination analysis outputs, Preseq curves, the Qualimap directory, samtools stats, the custom coverage and reads processing summary tables, the mapDamage2 output directory, and the DeDup JSON summary. Each input type is only requested if its corresponding config toggle is enabled — disabled steps are silently omitted, keeping the report self-consistent.

---

## Module 3 — Dynamics Processing

This module quantifies the relative abundance and activity of transposable elements (TEs) and other genomic features across individuals. It works by mapping reads to a purpose-built reference consisting of TE sequences combined with single-copy genes (SCGs). The SCGs serve as a stable normalisation reference, making TE abundance estimates comparable across individuals regardless of differences in sequencing depth.

### Library Preparation

Two reference libraries are required: a **feature library** containing the TE (or other feature) sequences of interest, placed under `{species}/raw/dynamics/feature_library/`, and an **SCG library** containing sequences from single-copy conserved genes, placed under `{species}/raw/dynamics/scg/`. Both support `.fna`, `.fasta`, and `.fa` formats.

Before mapping, sequence headers in both libraries are standardised and suffixed to make TE and SCG sequences distinguishable in the alignments — `_fle` is appended to every feature library header, and `_scg` to every SCG header. The two processed libraries are then concatenated into a single combined FASTA, which is indexed with **BWA-MEM2** in preparation for mapping. Multiple feature libraries per species are supported; each produces an entirely independent set of downstream results.

### Mapping to the Combined Library

Merged per-individual reads (from Module 1) are mapped to the combined SCG + TE library using **BWA-MEM2**. After conversion from SAM to BAM, unmapped reads are immediately discarded to reduce file sizes. The filtered BAM is then sorted and indexed. Intermediate files are cleaned up automatically.

### Normalisation

For each individual, the pipeline calculates read coverage across both the SCG and TE sequences in the combined library. The SCG coverage serves as the normalisation factor, accounting for differences in sequencing depth and mapping efficiency between individuals. This normalisation is critical for making meaningful comparisons of TE abundance across a population. Per-individual coverage files are combined into a species-level table, and an R-based plot is generated showing normalised TE abundance across all individuals. The order and labels of individuals in the plot are derived from the individual list for the species.

### TEplotter Analysis

TEplotter provides a complementary view focusing on **site occupancy (SO)** — a measure of how consistently a given TE insertion site is represented across individuals. The analysis runs through a multi-step pipeline:

First, the sorted BAM is converted into a site-occupancy profile, using the combined library FASTA to determine sequence lengths. The raw SO values are then normalised, and occupancy statistics are estimated. The normalised profiles are converted into a plotable directory structure, processing all sequences and labelling outputs with the individual's identifier.

From the plotable directories, per-individual TE occupancy plots are generated. A second pass combines all individual plotable directories to produce a faceted species-level comparison plot, allowing side-by-side visual comparison of TE dynamics across all individuals simultaneously.

---

## Module 4 — Processing Summary

The summary module consolidates outputs from all three processing modules into cohesive **MultiQC** HTML reports at two levels: a per-individual report covering all QC and analytics for a single individual across all references, and a per-species overall report aggregating all individuals.

Custom data preparation scripts transform pipeline outputs — coverage breadth, coverage depth, and reads processing statistics — into TSV files formatted for MultiQC's custom content sections. Reads processing summaries are produced both as absolute values (raw, trimmed, quality-filtered, endogenous, deduplicated counts) and as a stacked representation (non-endogenous, duplicates, retained endogenous reads) suited to stacked bar chart visualisation in MultiQC. Both formats are additionally combined into species-level files for the overall report.

Because Qualimap and mapDamage2 output directories need to be co-located with other MultiQC inputs, they are copied into the summary folder structure before report generation. All inputs to MultiQC are conditionally requested based on which pipeline steps were enabled, so the reports always accurately reflect the steps that were actually run.

---

## Global Configuration and Pipeline Behaviour

The entire Past Forward aDNA Pipeline is controlled through a single `config/config.yaml` file. Species to be processed are defined as a top-level mapping, each with an optional display name, and all modules run for every species listed.

Each major processing step can be independently enabled or disabled. The key toggles and their defaults are:

- **`pipeline > raw_reads_processing > adapter_removal > execute`** — enable adapter removal (default `true`)
- **`pipeline > raw_reads_processing > quality_filtering > execute`** — enable quality filtering (default `true`)
- **`pipeline > reads_processing > quality_checking_merged > execute`** — enable FastQC on merged reads (default `true`)
- **`pipeline > raw_reads_processing > contamination_analysis > execute`** — enable contamination analysis (default `true`)
- **`pipeline > reference_processing > execute`** — enable the entire reference processing module (default **`false`**, must be explicitly set to `true`)
- **`pipeline > reference_processing > deduplication > execute`** — enable deduplication (default `true`)
- **`pipeline > reference_processing > damage_rescaling > execute`** — enable mapDamage2 rescaling (default `true`)
- **`pipeline > reference_processing > coverage_analysis > execute`** — enable coverage and BAM statistics (default `true`)

A global **`pipeline > global > skip_existing_files`** option (default `true`) checks for already-completed outputs before the pipeline starts and excludes them from the target list, preventing redundant reprocessing when resuming an interrupted run.

On every execution the pipeline logs extensive provenance information: timestamp, platform and OS details, Python and Snakemake versions, the active conda environment, the git commit hash of the pipeline code, the full command line used, all config file paths, and the complete loaded configuration. A minimum Snakemake version of 9.9.0 is enforced at startup.