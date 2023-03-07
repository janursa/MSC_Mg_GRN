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

from scripts.imports import GRN_DIR, MODELSELECTION_DIR, F_DE_protiens
from scripts.utils.links import plot_mean_weights
from scripts.utils import make_title_pretty

if __name__ == '__main__':
    studies = ['ctr', 'mg']
    methods = ['RF', 'ridge', 'portia']
    selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']

    nrows = 1
    ncols = 2
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(ncols * 3.2, nrows * 2.5))
    model_i = 0
    #- plot selected models
    for DE_type, _ in F_DE_protiens().items():
        for method in methods:
            model_name = '_'.join([DE_type, method])
            if model_name not in selected_models:  # only selected models
                continue
            links_stack = [pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=[0]) for study in studies]
            # print(links_stack[0].sort_values('Weight', ascending=False),'\n')
            if axes.shape[0] == 1:
                ax = axes[model_i]
            else:
                i = int(model_i/ncols)
                j = model_i%ncols
                ax = axes[i][j]
            #- calculate what weights are in terms of top quantile

            plot_mean_weights(ax, links_stack, make_title_pretty(model_name), studies)
            model_i+=1
    fig.savefig(os.path.join(GRN_DIR, f'mean_weights.pdf'))
    fig.savefig(os.path.join(GRN_DIR, f'mean_weights.png'), dpi=300, transparent=True)

