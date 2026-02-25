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
- Before running the pipeline, ensure you have an environment with Snakemake and the required dependencies installed.
- Required dependencies for pipeline processing will be installed automatically, except for the contamination analysis tools. Those have to be installed separately and their details need to be added to the config file.
  - The pipeline supports ECMSD for contamination analysis. Ensure ECMSD is configured in the `config.yaml` file under `contamination_analysis`.
  - The pipeline supports Centrifuge for contamination analysis. Ensure Centrifuge is configured in the `config.yaml` file under `contamination_analysis`.
- You need to add species details to the pipeline (config and files).
- your reads should be renamed according to the naming convention specifie dbelow

## Folder Structure

### Species Folders

The project contains folders for different species, which contain the raw data, processed data, and results for each species.

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

General reads processing data will be in either `processed`or `results`. 

Everything related to a reference will have a `<reference>` folder under `processed`or `results`. Typically, only the `results` folder will contain information required for further analyis. In case more information is required, the original files can often be found in the `processed` folder. 

Some exemptions include `*.sam` and unsorted `*.bam` files. These are deleted to save storrage space. Most other files are kept in order to allow reprocessing the pipeline from different points in case something fails. If a step should be repeated, the relevant files need to be deleted manually. 

#### RAW Reads Filenames

The pipeline expects input read files to follow a standardized naming convention:

```
<Individual>_[<FreeText>_]R<1/2>[_<FreeText>].fastq.gz
```

Following this convention ensures proper organization and automated processing within the pipeline.  

##### Filename Components:
- **`<Individual>`** – A unique identifier for the sample or individual.  
- **`<FreeText>`** – Any additional text or identifier that can be included in the filename. Typically, this is used to differentiate between different samples within the same individual. E.g. the original sample name.  
- **`R<1/2>`** – Indicates the read pair number, typically `R1` for the first read and `R2` for the second read. 
- **`.fastq.gz`** – The expected file extension, indicating compressed FASTQ format. Only `.fastq.gz` files are supported.

Notes:
 - This name must contain `_R1` and, if paired-end, `_R2`.
 - For paired-end data, the name of the reads must be identical except for `_R1` and `_R2`.
 - Individual names/IDs will be used to name the output files as well as in the reports and plots.

#### Example:
```
Bger1_326862_S37_R1_001.fastq.gz
```

# Configuration File Structure for aDNA Pipeline (`config.yaml`)

The `config.yaml` file is used to configure the aDNA pipeline. It contains settings such as project name, the species list and the pipeline stages and their process steps.

## Global Settings

* **project\_name**: Name of the aDNA project.

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

### Example `config.yaml`

```yaml
# config.yaml - Configuration file for aDNA pipeline
# This file contains settings for various stages of the pipeline

project_name: "aDNA_Project"

# Pipeline stages and their configurations
pipeline:

  # Global settings
  global:
    # When true, existing output files will be skipped to avoid re-computation (Default: true)
    skip_existing_files: true

  # Stages of the pipeline

  # Raw reads processing
  # Includes quality checking, adapter removal, quality filtering, merging, 
  # contamination analysis, and statistical analysis
  raw_reads_processing:
    # When true, this stage will be executed. (Default: true)
    execute: true

    # Sub-stages with their respective settings
    # Quality checking of raw reads
    quality_checking_raw:
      # When true, this sub-stage will be executed (Default: true)
      execute: true
    
    # Adapter removal from raw reads
    adapter_removal:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

      # Settings for adapter removal
      settings: 
        # Minimum quality score for adapter removal
        min_quality: 0
        # Minimum length of reads after adapter removal
        min_length: 0
        # Optional: Adapter sequences for read 1 and read 2
        # If not provided, fastp will try to identify adapters automatically
        adapters_sequences:
          # Adapter sequence for read 1
          r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
          # Adapter sequence for read 2
          r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
    
    # Quality checking of trimmed reads
    quality_checking_trimmed:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Quality filtering of trimmed reads
    quality_filtering:
      # When true, this sub-stage will be executed (Default: true)
      execute: true
      settings:
        # Minimum quality score for quality filtering
        min_quality: 15
        # Minimum length of reads after quality filtering
        min_length: 30

    # Quality checking of quality-filtered reads
    quality_checking_quality_filtered:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Quality checking of merged reads
    quality_checking_merged:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Contamination analysis
    contamination_analysis:
      # When true, this sub-stage will be executed (Default: true)
      execute: true
      tools:
        # ECMSD tool settings for contamination analysis
        ecmsd:
          # When true, this tool will be executed (Default: true)
          execute: true
          settings:
            # Optional: Path to the conda environment for ECMSD
            # If not provided, the default environment will be used
            #conda_env: "../../../../envs/ecmsd.yaml"
            # Path to the ECMSD executable
            # Curretnly, ecmsd can not be installed via conda. Provide the path to the shell script to run ECMSD.
            executable: "/path/to/ecmsd/shell/ECMSD.sh"
        # Centrifuge tool settings for contamination analysis
        centrifuge:
          # When true, this tool will be executed (Default: true)
          execute: true
          settings:
            # Optional: Path to the conda environment for Centrifuge
            # If not provided, the default environment will be used
            #conda_env: "../../../../envs/centrifuge.yaml"
            # Path to the Centrifuge index
            index: "/path/to/centrifuge_index"
    
    # Statistical analysis
    statistical_analysis:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

  # Reference processing
  reference_processing:
    # When true, this stage will be executed (Default: true)
    execute: true
    
    # Deduplication settings
    deduplication:
      # When true, this sub-stage will be executed (Default: true)
      execute: false
      settings:
        # To increase performance, deduplication will be done per cluster of contigs
        # Below settings define how the contigs will be clustered
        # Optional: Maximum number of contigs per cluster (Default: 10 if not specified)
        max_contigs_per_cluster: 10
        # Optional: Maximum number of contigs per cluster (Default: 500 if not specified)
        max_contigs_per_cluster: 500
    
    # Damage rescaling settings for mapDamage2
    damage_rescaling:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Damage analysis settings for mapDamage2
    damage_analysis:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Endogenous reads analysis settings
    endogenous_reads_analysis: 
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Coverage analysis settings
    coverage_analysis: 
      # When true, this sub-stage will be executed (Default: true)
      execute: true

# Species details
species:
  Bger:
    name: "Blatella germanica"
  Dsim:
    name: "Drosophila simulans"