"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, F_DE_proteins
from scripts.utils.links import normalize_links

def retreive_grn(DE_type, study, methods):
    links_stack = [] #all methods
    for method in methods:
        links_stack.append(pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'),index_col=[0]))
    return links_stack
if __name__ == '__main__':
    studies = ['ctr', 'mg']
    methods = ['RF', 'ridge', 'portia']
    selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    ensemble_dir = os.path.join(GRN_DIR, 'ensemble')
    if not os.path.isdir(ensemble_dir):
        os.makedirs(ensemble_dir)

    for study in studies:
        links_stack = [] # different models
        for idx, (DE_type, DE_proteins) in enumerate(F_DE_proteins().items()):
            for method in methods:
                model_name = '_'.join([DE_type, method])
                if model_name not in selected_models:  # only selected models
                    continue
                links_stack.append(
                    pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv')))

        for links in links_stack:
            print(links)
            print('\n')
        aa
        links_stack_n = [normalize_links(links) for links in links_stack]
        links_ensemble = links_stack_n[0]

        links_ensemble['Weight'] = np.mean([links['Weight'] for links in links_stack_n], axis=0)
        links_ensemble.to_csv(os.path.join(ensemble_dir, f'links_{study}.csv'), index=False)
