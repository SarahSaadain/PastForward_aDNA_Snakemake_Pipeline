# PastForward - An aDNA Snakemake Pipeline

This project contains a pipeline to analyze raw ancient data, obtained from the sequencing facility. The pipeline includes various Snakemake workflows to process, analyze, and generate reports on the sequence quality, which helps decide if an aDNA extraction and sequencing was successful, and further polishes the data for downstream analyses.

Note: This pipeline is still in the final stages of development. It can already be used for analysis but might still be subject to changes.

## Workflow Overview

Below is a simplified overview of the processing steps of the pipeline:

![Pipeline Overview](docs/img/aDNA_pipeline_process.svg)

For detailed information about the processing steps, see the [Process Overview](docs/process_overview.md) page.

## Setup Overview

The PastForward pipeline is implemented using Snakemake. Information about the setup as well as configuration options can be found in the [Setup Instructions](config/README.md).

For more information on Snakemake, see the [Snakemake website](https://snakemake.github.io).

## Running the Pipeline

The PastForward pipeline is implemented using Snakemake, a workflow management system. Snakemake ensures reproducibility and efficient execution of the pipeline.

### Running the Pipeline

To run the pipeline, navigate to the root directory containing the `workflow/Snakefile` and execute:

```bash
# minimum command to run the pipeline
#snakemake --cores <number_of_threads> --use-conda

# suggested command to run the pipeline
snakemake --cores <number_of_threads> --use-conda --keep-going –rerun-trigger mtime
```

Replace `<number_of_threads>` with the number of CPU threads you want to allocate for the pipeline.

**Note:** 
* The `--use-conda` flag enables the use of conda environments specified in the `Snakefile`.
* The `--keep-going` flag allows the pipeline to continue even if a rule fails. Somtimes the analysis of ECMSD fails due to issues with the input data (e.g. low quality reads or low coverage). In this case, the rest of the pipeline can still be executed.
* The number of threads can be adjusted using the `--cores` option when running Snakemake.
* The `–rerun-trigger mtime` flag ensures that the pipeline only re-runs rules if the input files have been modified since the last run.

Other useful flags:
* `--dryrun` or `-n` to simulate the execution of the pipeline without actually running it
* `--configfile <path_to_config.yaml>` to specify a custom config file
* `--rerun-incomplete` to re-run rules that failed or were cancelled in the previous run
* `--rerun-trigger` to specify which triggers to use for rerunning rules
  * Possible choices: code, input, mtime, params, software-env
  * Define what triggers the rerunning of a job. By default, all triggers are used, which guarantees that results are consistent with the workflow code and configuration. If you rather prefer the traditional way of just considering file modification dates, use `–rerun-trigger mtime`.
* `--touch` to touch output files (mark them up to date without really changing them) instead of running their commands. This is used to pretend that the rules were executed, in order to fool future invocations of snakemake. Note that this will only touch files that would otherwise be recreated by Snakemake (e.g. because their input files are newer). For enforcing a touch, combine this with –force, –forceall, or –forcerun. Note however that you lose the provenance information when the files have been created in reality. Hence, this should be used only as a last resort.

For more information on Snakemake command-line options, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

### Running the Pipeline in the Background

Depending on the size of the data, it may take some time to complete the pipeline. Thus it is recommended to run the pipeline in the background. You can do this by running the following command:

```bash
nohup snakemake --cores 40 --use-conda --keep-going --rerun-trigger mtime > pipeline.log 2>&1 &
```

### Restarting the Pipeline

Snakemake automatically tracks the state of the pipeline and will only re-run steps that are incomplete or outdated. If you want to restart the pipeline from the beginning, you can delete the relevant output files and re-run the pipeline.

If you want to restart the pipeline, because it has crashed or was terminated, you might need to use the `--rerun-incomplete` flag. This will re-run all incomplete steps, even if they have not been modified since the last run.

## Reports

The PastForward pipeline generates several MultiQC reports to provide a comprehensive summary of the quality control and analysis results at various stages of the workflow. These reports are essential for assessing the quality of the sequencing data and the results of the processing pipeline.

By leveraging the AI functionality in the MultiQC reports, you can use ai to interpret the results of the pipeline and make informed decisions about the quality of the data and the results of the analysis.

Note: Currently, the report functionality of snakemake is not available for the PastForward pipeline.

### BAM File MultiQC Reports

The BAM file MultiQC reports are a key output of the PastForward pipeline and are essential for downstream analysis and quality checking. These reports are generated at two levels:

1. **Reference-Level BAM File MultiQC Report**:
   - **Location**: `{species}/results/{reference}/analytics/{individual}_{reference}_multiqc.html`
   - **Description**: Summarizes the quality metrics of BAM files, including results from reads processing, contamination, coverage analysis, deduplication, and damage rescaling and additiional statistics. Provides a detailed report for each individual and reference.

2. **Individual-Level BAM File MultiQC Report**:
   - **Location**: `{species}/results/summary/{individual}_multiqc.html`
   - **Description**: Summarizes the quality metrics of BAM files, including results from reads processing, contamination, coverage analysis, deduplication, and damage rescaling and additiional statistics. Provides a detailed report for each individual across all references.

3. **Species-Level BAM File MultiQC Report**:
   - **Location**: `{species}/results/summary/{species}_multiqc.overall.html`
   - **Description**: Summarizes the quality metrics of BAM files, including results from reads processing, contamination, coverage analysis, deduplication, and damage rescaling and additiional statistics. Provides a detailed report for all individuals across all references.

### Additional MultiQC Reports for Reads

In addition to the BAM file reports, the pipeline generates MultiQC reports for raw, trimmed, quality-filtered, and merged reads. These reports provide more detailed insights into the quality of the sequencing reads and can be used for additional analysis if required. The locations of these reports are as follows:

1. **Raw Reads MultiQC Report**:
   - **Location**: `{species}/results/reads/{species}_multiqc_raw.html`
   - **Description**: Summarizes the quality metrics of raw sequencing reads.

2. **Trimmed Reads MultiQC Report**:
   - **Location**: `{species}/results/reads/{species}_multiqc_trimmed.html`
   - **Description**: Provides quality metrics for reads after adapter trimming.

3. **Quality-Filtered Reads MultiQC Report**:
   - **Location**: `{species}/results/reads/{species}_multiqc_quality_filtered.html`
   - **Description**: Details the quality metrics of reads after quality filtering.

4. **Merged Reads MultiQC Report**:
   - **Location**: `{species}/results/reads/{species}_multiqc_merged.html`
   - **Description**: Contains quality metrics for reads after merging paired-end reads.