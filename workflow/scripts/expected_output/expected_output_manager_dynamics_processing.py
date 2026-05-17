

import logging

def get_expected_output_dynamics_processing(species):
    all_inputs = []

    if config.get("pipeline", {}).get("dynamics", {}).get("execute", True) == False:
        logging.info(f"Skipping dynamics processing for {species}. Disabled in config.")
        return []
    
    try:
        feature_libraries = get_feature_library_ids_for_species(species)
    except Exception as e:
        logging.warning(f"No feature libraries defined for {species}. Skipping dynamics processing.")
        return []
    
    try:
        scgs = get_scg_library_ids_for_species(species)
    except Exception as e:
        logging.warning(f"No SCG libraries defined for {species}. Skipping dynamics processing.")
        return []
    
    # only support 1 scg library
    if len(scgs) > 1:
        raise ValueError(f"Multiple SCG libraries provided for {species}. Only one SCG library is supported.")

    individuals = get_individuals_for_species(species)

    all_inputs.append(f"{species}/results/dynamics/{species}_seqvista_stats_comparison.tsv")

    seqvista_settings = config.get("pipeline", {}).get("dynamics", {}).get("seqvista", {}).get("settings", {})
    individual_plots_mode = seqvista_settings.get("individual_plots", "plot")
    comparison_plots_mode = seqvista_settings.get("comparison_plots", "plot")

    for feature_library in feature_libraries:

        if config.get("pipeline", {}).get("dynamics", {}).get("seqvista", {}).get("execute", True) == True:
            #all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/{species}_estimation.combined.tsv")

            # Species-level stats (always produced when seqvista is enabled)
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/species_level/{species}_{feature_library}_stats_comparison.tsv")
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/species_level/{species}_{feature_library}_flagged_seqids.tsv")

            # Species-level comparison plots
            if comparison_plots_mode == "plot":
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/species_level/{species}_plots_facet/")
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/species_level/{species}_plotables_facet.tar.gz")
            elif comparison_plots_mode == "plotable_only":
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/species_level/{species}_plotables_facet.tar.gz")

            for individual in individuals:
                # Coverage files (always produced when seqvista is enabled)
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/individual_level/{individual}_coverage.tsv.gz")
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/individual_level/{individual}_coverage.normalized.tsv.gz")
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/individual_level/{individual}_coverage.normalized.stats.tsv")

                # Individual-level plots
                if individual_plots_mode == "plot":
                    all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/individual_level/{individual}_plotable.tar.gz")
                    all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/individual_level/{individual}_plots/")
                elif individual_plots_mode == "plotable_only":
                    all_inputs.append(f"{species}/results/dynamics/{feature_library}/seqvista/individual_level/{individual}_plotable.tar.gz")
    
        if config.get("pipeline", {}).get("dynamics", {}).get("pf_normalization", {}).get("execute", True) == True:
            
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/normalization/plots/")
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/normalization/{species}_normalized_coverage.combined.tsv")

    return all_inputs