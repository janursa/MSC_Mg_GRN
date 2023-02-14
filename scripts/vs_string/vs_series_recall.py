"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import ENRICH_DIR, VS_STRING_DIR, F_DE_data, GRN_DIR, F_DE_protiens, CALIBRATION_DIR
from scripts.utils.links import plot_match_counts_series, compare_network_string, format_links_string
from scripts.utils import calibration, serif_font, make_title_pretty

from geneRNI.evaluation import precision_recall_curve

if __name__ == '__main__':
    if not os.path.isdir(VS_STRING_DIR):
        os.makedirs(VS_STRING_DIR)
    n_repeat = 10

    DE_types = F_DE_data().keys()
    study='ctr'
    methods = ['RF', 'ridge', 'portia', 'arbitrary']

    ncols = 2
    nrows = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))

    links_string_dict = {}
    for DE_type in DE_types:
        links_string = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        links_string_dict[DE_type] = links_string

    serif_font()
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey', 'lightcoral']
    linestyles = ['-', '--', '-.', ':', '-', '--']

    for idx, DE_type in enumerate(DE_types):
        i = int(idx/ncols)
        j = idx%ncols
        ax = axes[i][j]
        golden_links = format_links_string(links_string_dict[DE_type], F_DE_protiens()[DE_type])
        protnames = F_DE_protiens()[DE_type]
        # - get the links, filter mask based on scores, and calculate recall
        for i, method in enumerate(methods):
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
            if method in ['ridge', 'RF']:
                best_scores, best_params = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
                protnames_left = np.asarray(protnames)[np.asarray(best_scores > 0)]
                filter_mask = links['Target'].isin(protnames_left)

            else:
                filter_mask = None
            _, recall, _ = precision_recall_curve(links, golden_links, filter_mask)

            label = method + f' (score = {str(round(np.mean(recall[:-1]), 2))})'
            ax.plot(np.linspace(0, 1, len(recall)), recall, label=label, color=colors[i], alpha=1, linewidth=2,
                        linestyle=linestyles[i])

        ax.set_title(DE_type)
        ax.legend(frameon=False)
        ax.set_ylabel('Recall')
        ax.set_xlabel('Thresholds')

    # fig_pr.savefig(os.path.join(VS_STRING_DIR, f'precision_recall.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(VS_STRING_DIR, f'recall.pdf'))




