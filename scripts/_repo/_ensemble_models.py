"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Callable


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))


from imports import MODELSELECTION_DIR, GRN_DIR
from utils.links import normalize_links

def retreive_links(models, studies):
    phase = models[0][0]
    links_stack_ctr_sample = [[] for i in range(len(studies))]
    for study_i, study in enumerate(studies):
        for model in models:
            imput, GRN_method = model[1], model[2]
            links = pd.read_csv(Path(GRN_DIR) / GRN_method / f'links_{phase}_{imput}_{study}.csv', index_col=False)
            links_stack_ctr_sample[study_i].append(links)
    return links_stack_ctr_sample

def ensemble_links(models, studies):
    links_stack_studies = retreive_links(models, studies)
    links_ensemble_studies = []
    links_template = links_stack_studies[0][0].copy()
    for study_i, study in enumerate(studies):
        links_stack_n = [normalize_links(links) for links in links_stack_studies[study_i]]
        links_template['Weight'] = np.mean([links['Weight'] for links in links_stack_n], axis=0)
        links_ensemble_studies.append(links_template)

    return links_ensemble_studies
def output_links(phase, studies, early_model_ensembles):
    for study, links_ensemble in zip(studies, early_model_ensembles):
        links_ensemble.to_csv(Path(MODELSELECTION_DIR) / f'links_ensemble_{phase}_{study}.csv', index_col=False)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    studies = args.studies

    #- load names of selected models
    shortlisted_model_names = np.loadtxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), dtype=str, delimiter=",")
    #- divide the models into early and late
    shortlisted_model = [name.split('_') for name in shortlisted_model_names]
    early_models = list(filter(lambda item: item[0]=='early', shortlisted_model))
    late_models = list(filter(lambda item: item[0] == 'late', shortlisted_model))
    #- obtain ensemble links for early and late, for both studies
    early_model_ensembles = ensemble_links(early_models, studies)
    late_model_ensembles = ensemble_links(late_models, studies)
    #- output the ensemble links for both studies
    output_links('early', studies, early_model_ensembles)
    output_links('late', studies, late_model_ensembles)