####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def combine_teplotters_for_species_input_coverage_files(wildcards):
    """
    Return the full path to the FASTQ file corresponding to this sample
    from the config.
    """
    species = wildcards.species
    feature_library = wildcards.feature_library

    individuals = get_individuals_for_species(species)

    list_of_teplotter_files_of_individuals = []

    for individual in individuals:
        list_of_teplotter_files_of_individuals.append(f"{species}/results/dynamics/{feature_library}/teplotter/{individual}_estimation.tsv")
    
    if not list_of_teplotter_files_of_individuals:
        raise ValueError(f"No teplotter files could be determined for species {species}.")

    return list_of_teplotter_files_of_individuals

####################################################
# Snakemake rules
####################################################

rule determine_teplotter_of_individual_bam_to_so:
    input:
        bam="{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam",
        fasta="{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
    output:
        coverage="{species}/results/dynamics/{feature_library}/teplotter/{individual}_coverage.tsv"
    log:
        "{species}/results/dynamics/{feature_library}/teplotter/{individual}_bam2so.log"
    conda:
        "../../envs/python.yaml"
    shell:
        """
        python workflow/scripts/dynamics/teplotter/bam2so.py --infile {input.bam} --fasta {input.fasta} --outfile {output.coverage} 2> {log}
        """

rule normalize_teplotter_of_individual:
    input:
        coverage="{species}/results/dynamics/{feature_library}/teplotter/{individual}_coverage.tsv"
    output:
        normalized="{species}/results/dynamics/{feature_library}/teplotter/{individual}_coverage.normalized.tsv"
    conda:
        "../../envs/python.yaml"
    shell:
        """
        python workflow/scripts/dynamics/teplotter/normalize-so.py --so {input.coverage} --outfile {output.normalized}
        """

rule estimate_teplotter_of_individual:
    input:
        coverage="{species}/results/dynamics/{feature_library}/teplotter/{individual}_coverage.tsv"
    output:
        estimation="{species}/results/dynamics/{feature_library}/teplotter/{individual}_estimation.tsv"
    conda:
        "../../envs/python.yaml"
    shell:
        """
        python workflow/scripts/dynamics/teplotter/estimate-so.py --so {input.coverage} --outfile {output.estimation}
        """

rule prepare_teplotter_visualization_of_individual:
    input:
        coverage="{species}/results/dynamics/{feature_library}/teplotter/{individual}_coverage.normalized.tsv",
    output:
        plotable=directory("{species}/results/dynamics/{feature_library}/teplotter/{individual}_plotable")
    conda:
        "../../envs/python.yaml"
    shell:
        """
        python workflow/scripts/dynamics/teplotter/so2plotable.py \
            --so {input.coverage} \
            --outdir {output.plotable} \
            --seq-ids ALL \
            --sample-id {wildcards.individual}
        """

rule run_teplotter_visualization_of_individual:
    input:
        "{species}/results/dynamics/{feature_library}/teplotter/{individual}_plotable"
    output:
        directory("{species}/results/dynamics/{feature_library}/teplotter/{individual}_plots")
    conda:
        "../../envs/python_and_r.yaml"
    params:
        log_threshhold = 25
    shell:
        """
        python workflow/scripts/dynamics/teplotter/run_plotable.py --folder {input} --outdir {output}  --log {params.log_threshhold}
        """

rule run_teplotter_visualization_of_species:
    input:
        lambda wildcards: expand(
            "{species}/results/dynamics/{feature_library}/teplotter/{individual}_plotable", 
            species=wildcards.species,
            feature_library=wildcards.feature_library,
            individual=get_individuals_for_species(wildcards.species))
    output:
        directory("{species}/results/dynamics/{feature_library}/teplotter/{species}_plots_facet")
    conda:
        "../../envs/python_and_r.yaml"
    params:
        log_threshhold = 25
    shell:
        """
        python workflow/scripts/dynamics/teplotter/run_plotable.py --folders {input} --outdir {output} --log {params.log_threshhold}
        """
