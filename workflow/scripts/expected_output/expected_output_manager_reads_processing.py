# =================================================================================================
#     Input Manager Utility Functions for pastForward Pipeline
# =================================================================================================
# This script provides helper functions for managing input files and sample metadata in the pipeline.
# Functions include file discovery, reference handling, and sample identification.

import os
import logging

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (raw reads)
def get_expected_output_fastqc_raw(species):

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_raw", {}).get("execute", True) == False:
        logging.info(f"Skipping FastQC for raw reads for {species}. Disabled in config.")
        return []

    files = get_r1_read_files_for_species(species)

    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','')
        all_inputs.append(f"{species}/results/reads/reads_raw/fastqc/{filename}_raw_fastqc.html")
        # Add R2 if exists
        if os.path.exists(raw_file.replace("_R1", "_R2")):
            all_inputs.append(f"{species}/results/reads/reads_raw/fastqc/{filename.replace('_R1', '_R2')}_raw_fastqc.html")
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (adapter trimmed reads)
def get_expected_output_fastqc_trimmed(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_trimmed", {}).get("execute", True) == False:
        logging.info(f"Skipping FastQC for trimmed reads for {species}. Disabled in config.")
        return []

    all_inputs = []
    for sample in get_sample_ids_for_species(species):
        all_inputs.append(f"{species}/results/reads/reads_trimmed/fastqc/{sample}_trimmed_fastqc.html")
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (adapter removed reads)
def get_expected_output_fastqc_quality_filtered(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_quality_filtered", {}).get("execute", True) == False:
        logging.info(f"Skipping FastQC for quality filtered reads for {species}. Disabled in config.")
        return []

    all_inputs = []
    for sample in get_sample_ids_for_species(species):
        all_inputs.append(f"{species}/results/reads/reads_quality_filtered/fastqc/{sample}_quality_filtered_fastqc.html")
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (merged reads)
def get_expected_output_fastqc_merged(species):
    all_inputs = []
    for individual in get_individuals_for_species(species):
        all_inputs.append(f"{species}/results/reads/reads_merged/fastqc/{individual}_merged_fastqc.html")
    return all_inputs

#-----------------------------------------------------------------------------------------------
# Get expected output file paths for contamination analysis (ECMSD)
def get_expected_output_contamination_ecmsd(species):  

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("contamination_analysis", {}).get("tools", {}).get("ecmsd", {}).get("execute", True) == False:
        logging.info(f"Skipping contamination analysis with ECMSD for {species}. Disabled in config.")
        return []

    expected_outputs = []

    for individual in get_individuals_for_species(species):
        for sample in get_samples_for_species_individual(species, individual):
            expected_outputs.append(f"{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary.txt")
    
    return expected_outputs

#-----------------------------------------------------------------------------------------------
# Get expected output file paths for contamination analysis (Centrifuge)
def get_expected_output_contamination_centrifuge(species):

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("contamination_analysis", {}).get("tools", {}).get("centrifuge", {}).get("execute", True) == False:
        logging.info(f"Skipping contamination analysis with Centrifuge for {species}. Disabled in config.")
        return []

    expected_outputs = []

    for individual in get_individuals_for_species(species):
        for sample in get_samples_for_species_individual(species, individual):
            expected_outputs.append(f"{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_report.tsv")
            expected_outputs.append(f"{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_taxon_counts.tsv")
            expected_outputs.append(f"{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_output.tsv.gz")

    return expected_outputs

#-----------------------------------------------------------------------------------------------
# Get expected output file paths for contamination analysis (all tools)
def get_expected_output_contamination(species):  

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("contamination_analysis", {}).get("execute", True) == False:
        logging.info(f"Skipping contamination analysis for {species}. Disabled in config.")
        return []

    expected_outputs = []

    # Add ECMSD contamination analysis outputs
    expected_outputs += get_expected_output_contamination_ecmsd(species)

    # Add Centrifuge contamination analysis outputs
    expected_outputs += get_expected_output_contamination_centrifuge(species)

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for MultiQC reports
def get_expected_output_multiqc(species):
    
    expected_outputs = []

    # Add MultiQC reports for different read processing stages
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_raw", {}).get("execute", True) == True:
        expected_outputs.append(f"{species}/results/reads/{species}_multiqc_raw.html")
    else:
        logging.info(f"Skipping MultiQC report for raw reads for {species}. Disabled in config.")

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_trimmed", {}).get("execute", True) == True:
        expected_outputs.append(f"{species}/results/reads/{species}_multiqc_trimmed.html")
    else:
        logging.info(f"Skipping MultiQC report for trimmed reads for {species}. Disabled in config.")

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_quality_filtered", {}).get("execute", True) == True:
        expected_outputs.append(f"{species}/results/reads/{species}_multiqc_quality_filtered.html")
    else:
        logging.info(f"Skipping MultiQC report for quality filtered reads for {species}. Disabled in config.")
    
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_merged", {}).get("execute", True) == True:
        expected_outputs.append(f"{species}/results/reads/{species}_multiqc_merged.html")
    else:
        logging.info(f"Skipping MultiQC report for merged reads for {species}. Disabled in config.")

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for read count plots
def get_expected_output_reads_plots(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("statistical_analysis", {}).get("execute", True) == False:
        logging.info(f"Skipping read plots generation for {species}. Disabled in config.")
        return []
    
    expected_outputs = []

    expected_outputs.append(f"{species}/results/reads/plots/{species}_read_counts.png")
    expected_outputs.append(f"{species}/results/reads/plots/{species}_read_counts_comparison_by_individual.png")

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for adapter removal
def get_expected_output_adapter_removal(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("adapter_removal", {}).get("execute", True) == False:
        logging.info(f"Skipping adapter removal for {species}. Disabled in config.")
        return []
    
    expected_outputs = []
    for sample in get_sample_ids_for_species(species):
        expected_outputs.append(f"{species}/processed/reads/reads_trimmed/{sample}_trimmed_final.fastq.gz")
        
    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for quality filtering
def get_expected_output_quality_filtering(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_filtering", {}).get("execute", True) == False:
        logging.info(f"Skipping quality filtering for {species}. Disabled in config.")
        return []
    
    expected_outputs = []
    for sample in get_sample_ids_for_species(species):
        expected_outputs.append(f"{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered_final.fastq.gz")
        
    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for read merging
def get_expected_output_read_merging(species):
    
    # read merging is should always be the result of raw read processing, so we don't check config for this step.
    expected_outputs = []
    for individual in get_individuals_for_species(species):
        expected_outputs.append(f"{species}/processed/reads/reads_merged/{individual}.fastq.gz")

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get all expected output file paths for raw read processing
def get_expexted_output_raw_read_processing(species):

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("execute", True) == False:
        logging.info(f"Skipping raw read processing for {species}. Disabled in config.")
        return []
    
    # Add FastQC outputs for raw reads
    expected_outputs = []

    expected_outputs += get_expected_output_adapter_removal(species)
    expected_outputs += get_expected_output_quality_filtering(species)
    expected_outputs += get_expected_output_read_merging(species)

    # Add MultiQC reports for different read processing stages
    expected_outputs += get_expected_output_multiqc(species)

    # Add ECMSD contamination analysis outputs
    expected_outputs += get_expected_output_contamination(species)

    # Add summary plots for read count comparisons
    expected_outputs += get_expected_output_reads_plots(species)

    return expected_outputs