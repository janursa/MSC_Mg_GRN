import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List
import matplotlib.pyplot as plt


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_model_name_2_method_and_DE_type, VSA_DIR, F_selected_models, \
    F_DE_genenames, F_protnames_to_genenames
from common_tools.role_analysis import role_analysis, determine_top_role_change, determine_critical_role_change
from VSA.vsa_aux import retreive_links_with_genenames
from geneRNI.evaluation import to_matrix

def check_validity(links_studies):
    for links, name in zip(links_studies, ['ctr', 'mg']):
        matrix = to_matrix(links)
        row = matrix.sum(axis=1)
        col = matrix.sum(axis=0)
        print(name)
        # print(row)
        # print(col)
        print(f'row -> std {row.std()}')
        print(f'col -> std {col.std()}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--top_quantile_role_change', type=float, default=.75,
                        help="Top quantile of biggest distance in the role change from ctr to sample")
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    top_quantile_role_change = args.top_quantile_role_change
    # - create VSA dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)
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
        top_role_change.to_csv(Path(VSA_DIR) / f'top_role_change_{model_name}.csv', index=False)
        critical_role_change.to_csv(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv', index=False)


