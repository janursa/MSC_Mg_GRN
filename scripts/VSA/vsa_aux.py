"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import json
import matplotlib.pyplot as plt
from typing import List, Dict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_model_name_2_method_and_DE_type, VSA_DIR, GRN_DIR, change_protnames_to_genenames_in_links
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data

def retreive_links_with_genenames(model_name, studies):
    """Reads links for selected model and both studies, convert the links annotation from protname
    to genename and returns links"""

    # - break the model's name into DE_type format
    method, DE_type = F_model_name_2_method_and_DE_type(model_name)
    period, imput = DE_type.split('_')

    # - VSA analysis for both studies
    links_studies = []
    for study in studies:
        links = read_from_file(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
        #- change the names from protnames to genenames
        links = change_protnames_to_genenames_in_links(links)
        links_studies.append(links)
    return links_studies

def output_file_name(noise_runs_dir,noise_type, model_name, std):
    return f'{noise_runs_dir}/results_{noise_type}_{model_name}_{std}.json'


def retrieve_VSA_results(model_name):
    """ Get the results of VSA analysis
    """
    top_role_change = read_from_file(Path(VSA_DIR) / f'top_role_change_{model_name}.csv')
    critical_role_change = read_from_file(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')
    target_genes = np.concatenate((top_role_change['Entry'].to_numpy(), critical_role_change['Entry'].to_numpy()))
    # target_genes = list(set(target_genes))
    return target_genes
