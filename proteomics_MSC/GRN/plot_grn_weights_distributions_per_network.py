"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import pandas as pd
import argparse
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from _imports import GRN_DIR, F_DE_protnames, F_selected_models, study_colors, PLOT_WEIGHTS_DIR, make_title_pretty
from common_tools.links import plot_grn_weights_distributions_per_network

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--GRN_methods', nargs='+', default=['RF', 'ridge', 'portia'])
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    methods = args.GRN_methods

    selected_models = F_selected_models()

    nrows = 1
    ncols = len(selected_models)
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(ncols * 3.2, nrows * 2.5))
    model_i = 0
    #- plot selected models
    for DE_type, _ in F_DE_protnames().items():
        for method in methods:
            model_name = '_'.join([DE_type, method])
            if model_name not in selected_models:  # only selected models
                continue
            links_stack = [pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=[0]) for study in studies]
            ax = axes[model_i]
            # - calculate what weights are in terms of top quantile
            plot_grn_weights_distributions_per_network(ax, links_stack, make_title_pretty(model_name), studies, study_colors)
            model_i+=1
    fig.savefig(os.path.join(PLOT_WEIGHTS_DIR, f'mean_weights.pdf'))
    fig.savefig(os.path.join(PLOT_WEIGHTS_DIR, f'mean_weights.png'), dpi=300, transparent=True)

