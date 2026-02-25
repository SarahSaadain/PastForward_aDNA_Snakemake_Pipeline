

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

    for feature_library in feature_libraries:

        if config.get("pipeline", {}).get("dynamics", {}).get("teplotter", {}).get("execute", True) == True:
            #all_inputs.append(f"{species}/results/dynamics/{feature_library}/teplotter/{species}_estimation.combined.tsv")
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/teplotter/{species}_plots_facet/")
            for individual in individuals:
                all_inputs.append(f"{species}/results/dynamics/{feature_library}/teplotter/{individual}_plots/")
    
        if config.get("pipeline", {}).get("dynamics", {}).get("pf_normalization", {}).get("execute", True) == True:
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/normalization/plots/")
            all_inputs.append(f"{species}/results/dynamics/{feature_library}/normalization/{species}_normalized_coverage.combined.tsv")

    return all_inputs