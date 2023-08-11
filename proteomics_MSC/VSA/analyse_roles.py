import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List
import matplotlib.pyplot as plt
import numpy as np

from proteomics_MSC.imports import VSA_DIR, VSA_VERBOSE_DIR
from proteomics_MSC.common_tools import F_model_name_2_method_and_DE_type, F_selected_models, F_DE_genenames, F_protnames_to_genenames
from proteomics_MSC.common_tools.role_analysis import role_analysis, determine_top_role_change, determine_critical_role_change
from vsa_aux import retreive_links_with_genenames
from geneRNI.evaluation import to_matrix

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)


if __name__ == '__main__':

    studies = config['studies']
    top_quantile_role_change = config['top_quantile_role_change']
    # - create VSA dir
    if not os.path.isdir(VSA_VERBOSE_DIR):
        os.makedirs(VSA_VERBOSE_DIR)
    # - load names of selected models
    selected_models = F_selected_models()

    # - the main function for VSA analysis.
    for model_name in selected_models:
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        genenames = F_DE_genenames()[DE_type]
        # - VSA analysis for both studies
        vsa_results_studies = []
        links_studies = retreive_links_with_genenames(model_name, studies)
        # - check the validity of links
        # check_validity(links_studies)
        # - analyse VSA
        for study_i, _ in enumerate(studies):
            links = links_studies[study_i]
            vsa_results = role_analysis(links, genenames)
            vsa_results_studies.append(vsa_results)
        # - determine top role change from ctr to sample
        top_role_change = determine_top_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                    top_quantile=top_quantile_role_change)
        # - determine those proteins with a critical role change from ctr to sample
        critical_role_change = determine_critical_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                              target_role=3)
        # - write VSA analysis for ctr and sample to files
        for study_i, study in enumerate(studies):
            vsa_results_studies[study_i].to_csv(Path(VSA_DIR) / f'vsa_{model_name}_{study}.csv', index=False)
        # - write top role changes and critical role changes to file
        to_save = Path(VSA_VERBOSE_DIR) / f'top_role_change_{model_name}.csv'
        top_role_change.to_csv(to_save, index=False)
        to_save = Path(VSA_VERBOSE_DIR) / f'critical_role_change_{model_name}.csv'
        critical_role_change.to_csv(to_save, index=False)

        # - aggregated genes names
        top_genes = np.array(top_role_change['Entry'].tolist() + critical_role_change['Entry'].tolist())
        to_save = f'{VSA_VERBOSE_DIR}/top_genes_{model_name}.csv'
        np.savetxt(to_save, top_genes, fmt='%s')
        print(f'output -> {to_save}')



