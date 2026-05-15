####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def clean_scg_library_name_input(wildcards):
    """
    Return the full path to the FASTA file corresponding to this sample
    from the config.
    """
    species = wildcards.species
    scg_library = wildcards.scg_library

    scg_library_path = get_scg_library_file_for_species_and_library(species, scg_library)

    if not scg_library_path:
        raise ValueError(f"No SCG library file could be determined for species {species} and library {scg_library}.")

    return scg_library_path

def clean_feature_library_name_input(wildcards):
    """
    Return the full path to the FASTA file corresponding to this sample
    from the config.
    """
    species = wildcards.species
    feature_library = wildcards.feature_library

    feature_library_path = get_feature_library_file_for_species_and_library(species, feature_library)

    if not feature_library_path:
        raise ValueError(f"No feature library file could be determined for species {species} and library {feature_library}.")

    return feature_library_path

####################################################
# Snakemake rules
####################################################
rule clean_feature_library_name:
    input:
        clean_feature_library_name_input
    output:
        temp("{species}/processed/dynamics/lib/feature_library/{feature_library}.clean.fasta")
    message: "Preparing TE library for {wildcards.species}"
    shell:
        """
        ln -sf $(realpath {input}) {output}
        """

rule prepare_feature_library:
    input:
        "{species}/processed/dynamics/lib/feature_library/{feature_library}.clean.fasta"
    output:
        temp("{species}/processed/dynamics/lib/feature_library/{feature_library}.suffixed.fasta")
    message: "Preparing TE library for {wildcards.species}"
    shell:
        # remove trailing whitespace from headers and append _fle to each header
        """
        sed -E '/^>/ s/[[:space:]]//g; /^>/ s/(_fle)?$/_fle/' {input} > {output}
        """

rule clean_scg_library_name:
    input:
        clean_scg_library_name_input
    output:
        temp("{species}/processed/dynamics/lib/scg/{scg_library}.clean.fasta")
    message: "Preparing SCG library for {wildcards.species}"
    shell:
        """
        ln -sf $(realpath {input}) {output}
        """

rule prepare_scg_library:
    input:
        "{species}/processed/dynamics/lib/scg/{scg_library}.clean.fasta"
    output:
        temp("{species}/processed/dynamics/lib/scg/{scg_library}.suffixed.fasta")
    message: "Preparing SCG library for {wildcards.species}"
    shell:
        # remove trailing whitespace from headers and append _scg to each header
        """
        sed -E '/^>/ s/[[:space:]]//g; /^>/ s/(_scg)?$/_scg/' {input} > {output}
        """
        
rule combine_scg_and_ref_library:
    input:
        scg= lambda wildcards: f"{wildcards.species}/processed/dynamics/lib/scg/{get_scg_library_ids_for_species(wildcards.species)[0]}.suffixed.fasta",
        fl="{species}/processed/dynamics/lib/feature_library/{feature_library}.suffixed.fasta"
    output:
        library="{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
    message: "Concatenating SCG and Feature libraries for {wildcards.species}"
    shell:
        """
        cat {input.scg} {input.fl} > {output.library}
        """