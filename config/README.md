# Setup Overview

## Install Snakemake
To install Snakemake, you can use conda, which is a package manager that simplifies the installation of software and its dependencies. You can create a new conda environment for Snakemake and install it using the following commands:

```bash
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
conda activate snakemake
snakemake --help
```

Refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more installation options and details.

## Setup Instructions
- Before running the pipeline, ensure you have an environment with Snakemake and it is activated.
- You need to add species details to the pipeline (config and files).
- Your reads should be renamed according to the naming convention specified below.

## Folder Structure

### Species Folders

The project contains folders for different species, which each contain the raw data, processed data, and results for the particular species.

The species folders should be placed in the root folder of your pipeline.

#### Providing Raw Data
The pipeline supports automatically moving the raw reads to the `<species>/raw/reads/` folder as well as the reference to the `<species>/raw/ref/` folder. Simply provide the files in the `<species>` folder. Alternatively, you can manually move the files to the respective folders.
  - provide the raw reads in `<species>/raw/reads/` folder
  - provide the reference in `<species>/raw/ref/` folder

When adding a new species, make sure to 
- the species folder should be placed in the root folder of your pipeline
- add the folder name should match the species key which is defined in `config.yaml` below `species:` 

#### Folder Structure

All other folders will be created and populated automatically

- Folder `<species>/processed/` contains the intermediary files during processing. Most of these files are marked as temporary and will be deleted at the end of the pipeline. Some files are kept to allow reprocessing the pipeline from different points in case something fails.
- Folder `<species>/results/` contains the final results and reports. 

Everything related to a reference will have a `<reference>` folder under `processed`or `results`. Typically, only the `results` folder will contain information required for further analyis. In case more information is required, the original files can often be found in the `processed` folder. 

Some exemptions include `*.sam` and unsorted `*.bam` files. These are deleted to save storage space. Most other files are kept in order to allow reprocessing the pipeline from different points in case something fails. If a step should be repeated, the relevant files need to be deleted manually. 

#### RAW Reads Filenames

The pipeline expects input read files to follow a standardized naming convention:

```
<Individual>_[<FreeText>_]R<1/2>[_<FreeText>].fastq.gz
```

Following this convention ensures proper organization and automated processing within the pipeline.  

##### Filename Components:
- **`<Individual>`** – A unique identifier for the sample or individual.  
- **`<FreeText>`** – Any additional text or identifier that can be included in the filename. Typically, this is used to differentiate between different samples within the same individual, e.g. the same sample was extracted twice using different protocols.
- **`R<1/2>`** – Indicates the read pair number, typically `R1` for the first read and `R2` for the second read. If the data is single-end, only `R1` should be present.
- **`.fastq.gz`** – The expected file extension, indicating compressed FASTQ format. Only `.fastq.gz` files are supported.

#### Example:
```
Dmel01_DabneyProtocol_R1_001.fastq.gz
```

# Configuration File Structure for aDNA Pipeline (`config.yaml`)

The `config.yaml` file is used to configure the aDNA pipeline. It contains settings such as project name, the species list and the pipeline stages and their process steps.

All pipeline stages are enabled by default, so a minimal config containing only the project name and species list is sufficient to run the pipeline without any further changes:

```yaml
project_name: "pastForward_Project"

species:
  Dmel:
    name: "Drosophila melanogaster"
```

To adjust any pipeline settings, open [config_designer.html](config_designer.html) in a browser. The interactive Config Designer guides you through all available options and exports a ready-to-use `config.yaml`.

## Global Settings

* **project\_name**: Name of the project.

## Pipeline Settings

Defines the overall pipeline behavior, including execution controls and process details.

### Pipeline Stages and Process Steps

* The pipeline is broken into **stages** (e.g., `raw_reads_processing`, `reference_processing`).
* Each stage contains multiple **process steps** (e.g., `adapter_removal`, `deduplication`, ...).
* Both stages and process steps can be controlled with `execute: true/false` flags to enable or disable them.
* Some process steps include additional configurable settings (e.g., adapter sequences, database paths, ...).
* If an enabled process step requires data from a previous stage which is disabled in the config, the pipeline will execute the disabled process step anyway.

### Important Defaults

* You **do not need to specify all stages or process steps** explicitly.
* Any **stage or process step not provided in the config defaults to `execute: true`** and will be executed.

### Global Pipeline Settings

| Setting | Default | Description |
|---|---|---|
| `pipeline.global.skip_existing_files` | `true` | When true, existing output files are skipped to avoid re-computation. |

### Stage: `raw_reads_processing`

Quality checking, adapter removal, quality filtering, merging, contamination analysis, and statistical analysis of raw reads.

#### `adapter_removal`

| Setting | Default | Description |
|---|---|---|
| `settings.min_quality` | `0` | Minimum base quality score for adapter trimming. |
| `settings.min_length` | `0` | Minimum read length after adapter removal. |
| `settings.adapters_sequences.r1` | auto-detect | Adapter sequence for read 1. If omitted, fastp detects adapters automatically. |
| `settings.adapters_sequences.r2` | auto-detect | Adapter sequence for read 2. If omitted, fastp detects adapters automatically. |
| `settings.extra_params` | — | Optional extra parameters passed directly to fastp. |

#### `quality_filtering`

| Setting | Default | Description |
|---|---|---|
| `settings.min_quality` | `15` | Minimum base quality score for quality filtering. |
| `settings.min_length` | `30` | Minimum read length after quality filtering. |

#### `contamination_analysis`

Both tools operate on quality-filtered reads and can be toggled independently.

**ECMSD** — maps reads against a curated mitochondrial reference database.

| Setting | Default | Description |
|---|---|---|
| `tools.ecmsd.settings.database` | — | Path to the ECMSD database folder. If omitted, the pipeline auto-creates a database at `resources/ecmsd_database` via `ECMSD --create-db`. |
| `tools.ecmsd.settings.Binsize` | `1000` | Reference genome bin size for coverage calculation. |
| `tools.ecmsd.settings.RMUS_threshold` | `0.15` | Minimum Relative Mapping Uniqueness Score for a taxon to be reported. |
| `tools.ecmsd.settings.mapping_quality` | `20` | Minimum mapping quality score to include a read. |
| `tools.ecmsd.settings.taxonomic_hierarchy` | `species` | Taxonomic level at which to aggregate and report results. Options: `species`, `genus`, `family`, `order`. |

**Centrifuge** — k-mer-based taxonomic classification against a user-provided database.

| Setting | Default | Description |
|---|---|---|
| `tools.centrifuge.settings.include_human_taxid` | `false` | When true, the human taxid is included in the Centrifuge analysis. |
| `tools.centrifuge.settings.index` | — | Optional path to the Centrifuge index prefix. If omitted, the default index will be downloaded automatically. |
| `tools.centrifuge.settings.conda_env` | — | Optional path to a custom conda environment for Centrifuge. |

### Stage: `reference_processing`

Mapping, deduplication, damage analysis, coverage — runs per individual per reference.

#### `mapping`

Merged per-individual reads are mapped to the reference genome. Mapping always runs when Reference Processing is enabled.

| Setting | Default | Description |
|---|---|---|
| `settings.mapper` | `bwa-mem2` | Mapper to use. Options: `bwa-aln` (classic seed-and-extend, recommended for short aDNA reads <70 bp), `bwa-mem2` (faster, for longer reads), `minimap2` (versatile, uses `-ax sr` preset for short reads). |
| `settings.mapper_extra_params` | — | Optional extra parameters passed directly to the mapper. For `bwa-aln`, defaults to `-n 0.01 -k 2 -l 1024 -o 2` (Oliva et al. 2021). |

#### `deduplication`

Removes PCR and sequencing duplicates using DeDup. Default: **off**.

| Setting | Default | Description |
|---|---|---|
| `settings.min_contigs_per_cluster` | `10` | Minimum number of contigs grouped into a cluster. Small contigs below this count are merged together before deduplication. |
| `settings.max_contigs_per_cluster` | `500` | Maximum number of contigs grouped per deduplication cluster. Lower values use less memory but increase runtime. Reduce (e.g. to 100) only for large, highly fragmented reference genomes. |

#### `filter_unmapped_reads`

Optionally removes or extracts reads that did not map to the reference. Default: **off**.

| Setting | Default | Description |
|---|---|---|
| `settings.action` | `remove` | What to do with unmapped reads: `remove` — write a mapped-reads-only BAM; `extract_fastq` — write unmapped reads to a compressed FASTQ; `extract_fasta` — write unmapped reads to a compressed FASTA. |

#### Other `reference_processing` steps

| Step | Default | Description |
|---|---|---|
| `damage_rescaling` | on | Profiles cytosine deamination and rescales base quality scores using mapDamage2. |
| `damage_analysis` | on | Visualises damage patterns and generates statistics for MultiQC reporting using mapDamage2. |
| `endogenous_reads_analysis` | on | Calculates the fraction of reads that mapped to the reference (endogenous content). |
| `coverage_analysis` | on | Computes coverage breadth and mean depth; runs Qualimap and Preseq for library complexity. |

### Stage: `dynamics`

TE and genomic feature abundance analysis — maps to a combined SCG + feature library for depth-normalised comparisons. Requires feature and SCG libraries placed in `{species}/raw/dynamics/`.

#### `mapping`

| Setting | Default | Description |
|---|---|---|
| `settings.mapper` | `bwa-mem2` | Mapper to use. Same options as `reference_processing.mapping`. |
| `settings.mapper_extra_params` | — | Optional extra parameters passed directly to the mapper. |

#### Other `dynamics` steps

| Step | Default | Description |
|---|---|---|
| `seqvista` | on | Generates per-position coverage profiles and per-individual TE occupancy plots. |
| `pf_normalization` | on | Normalises TE abundance by SCG coverage for cross-individual comparability. |

### Stage: `summary_processing`

Consolidates all QC outputs into MultiQC HTML reports — one per individual and one per species.

### Example `config.yaml`

```yaml
# config.yaml - Configuration file for pastForward
# This file contains settings for various stages of the pipeline

project_name: "pastForward_Project"

# Pipeline stages and their configurations
pipeline:

  # Global settings
  global:
    # When true, existing output files will be skipped to avoid re-computation (Default: true)
    skip_existing_files: true

  # Raw reads processing
  raw_reads_processing:
    execute: true

    quality_checking_raw:
      execute: true

    adapter_removal:
      execute: true
      settings:
        min_quality: 0
        min_length: 0
        # Optional: custom adapter sequences (auto-detected by fastp if not provided)
        #adapters_sequences:
        #  r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        #  r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        # Optional: extra parameters for fastp
        #extra_params: ""

    quality_checking_trimmed:
      execute: true

    quality_filtering:
      execute: true
      settings:
        min_quality: 15
        min_length: 30

    quality_checking_quality_filtered:
      execute: true

    quality_checking_merged:
      execute: true

    contamination_analysis:
      execute: true
      tools:
        ecmsd:
          execute: true
          settings:
            # Optional: path to ECMSD database folder (auto-created if not set)
            #database: "resources/ecmsd_database"
            # Bin size for coverage calculation (Default: 1000)
            Binsize: 1000
            # Minimum RMUS for a taxon to be reported (Default: 0.15)
            RMUS_threshold: 0.15
            # Minimum mapping quality to include a read (Default: 20)
            mapping_quality: 20
            # Taxonomic level for results: "species" (default), "genus", "family", "order"
            taxonomic_hierarchy: "species"
        centrifuge:
          execute: true
          settings:
            # When true, human will be included in the Centrifuge analysis (Default: false)
            include_human_taxid: false
            # Optional: path to the Centrifuge index (downloaded automatically if not provided)
            #index: "/path/to/centrifuge_index"
            # Optional: custom conda environment for Centrifuge
            #conda_env: "../../../../envs/centrifuge.yaml"

    statistical_analysis:
      execute: true

  # Reference processing
  reference_processing:
    execute: true

    mapping:
      settings:
        # Options: "bwa-mem2" (default), "bwa-aln", "minimap2"
        mapper: "bwa-mem2"
        # Optional: extra parameters passed directly to the mapper
        #mapper_extra_params: ""

    deduplication:
      # Default: true — set to false only if library complexity is very low
      execute: true
      settings:
        # Minimum contigs per deduplication cluster (Default: 10)
        min_contigs_per_cluster: 10
        # Maximum contigs per deduplication cluster (Default: 500)
        max_contigs_per_cluster: 500

    damage_rescaling:
      execute: true

    filter_unmapped_reads:
      # Default: false — enable to remove or extract unmapped reads
      execute: false
      settings:
        # Options: "remove", "extract_fastq", "extract_fasta"
        action: "remove"

    damage_analysis:
      execute: true

    endogenous_reads_analysis:
      execute: true

    coverage_analysis:
      execute: true

  dynamics:
    execute: true

    mapping:
      settings:
        # Options: "bwa-mem2" (default), "bwa-aln", "minimap2"
        mapper: "bwa-mem2"
        # Optional: extra parameters passed directly to the mapper
        #mapper_extra_params: ""

    seqvista:
      execute: true

    pf_normalization:
      execute: true

  summary_processing:
    execute: true

# Species details
species:
  Dmel:
    name: "Drosophila melanogaster"
```
