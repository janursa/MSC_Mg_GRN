"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import numpy as np
import json
from typing import Dict, List

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_model_name_2_method_and_DE_type, GRN_DIR, VSA_DIR, F_DE_protnames, F_selected_models, F_DE_genenames, F_protnames_to_genenames
from utils.VSA import role_analysis, RolePlot, determine_top_role_change, determine_critical_role_change
from VSA.vsa_aux import retreive_links_with_genenames, plot_roles



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--top_quantile_role_change', type=float, default=.9, help="Top quantile of biggest distance in the role change from ctr to sample")
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    top_quantile_role_change = args.top_quantile_role_change
    #- create VSA dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)
    # - load names of selected models
    selected_models = F_selected_models()

    #- read the links and store it for selected models
    links_all: Dict[str, List[pd.DataFrame]] = retreive_links_with_genenames(selected_models, F_protnames_to_genenames(), studies)

    #- the main function for VSA analysis.
    for model_name in selected_models:
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        genenames = F_DE_genenames()[DE_type]
        #- VSA analysis for both studies
        vsa_results_studies = []
        links_studies = links_all[model_name]
        for study_i, _ in enumerate(studies):
            vsa_results = role_analysis(links_studies[study_i], genenames)
            vsa_results_studies.append(vsa_results)
        #- determine top role change from ctr to sample
        top_role_change = determine_top_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                       top_quantile=top_quantile_role_change)
        #- determine those proteins with a critical role change from ctr to sample
        critical_role_change = determine_critical_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                                     target_role=3)
        #- write VSA analysis for ctr and sample to files
        for study_i, study in enumerate(studies):
            vsa_results_studies[study_i].to_csv(Path(VSA_DIR) / f'vsa_{model_name}_{study}.csv', index=False)
        #- write top role changes and critical role changes to file
        top_role_change.to_csv(Path(VSA_DIR)/ f'top_role_change_{model_name}.csv', index=False)
        critical_role_change.to_csv(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv', index=False)

        #- plot roles for ctr and sample, top role changes and critical role changes in one plot
        fig = plot_roles(studies, vsa_results_studies[0], vsa_results_studies[1], top_role_change, critical_role_change, genenames)
        fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.pdf'))

