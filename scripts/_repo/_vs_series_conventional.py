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

from scripts.imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, F_DE_protiens
from scripts.utils.links import plot_match_counts_series, compare_network_string, compare_network_string_batch, create_random_links, format_links_string
from scripts.utils import serif_font, make_title_pretty

from geneRNI.evaluation import precision_recall_curve, calculate_auc_roc, PR_curve_gene, roc_curve

if __name__ == '__main__':
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)
    n_repeat = 10

    DE_types = F_DE_data().keys()
    study='ctr'
    methods = ['RF', 'ridge', 'portia']
    methods = ['ridge', 'portia'] #TODO: remove this

    ncols = 1
    nrows = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))
    fig_pr, axes_pr = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))
    fig_roc, axes_roc = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))

    links_string_dict = {}
    for DE_type in DE_types:
        links_string = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        links_string_dict[DE_type] = links_string

    serif_font()
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey', 'lightcoral']
    linestyles = ['-', '--', '-.', ':', '-', '--']
    for idx, method in enumerate(methods):
        ax = axes[idx]
        ax_pr = axes_pr[idx]
        ax_roc = axes_roc[idx]
        # - get the links
        for i, DE_type in enumerate(DE_types):
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)

            golden_links = format_links_string(links_string_dict[DE_type], F_DE_protiens()[DE_type])
            precision, recall, thresholds = precision_recall_curve(links, golden_links)
            fpr, tpr, _ = roc_curve(links, golden_links)
            # F1_score = np.sum(2 * (precision * recall) / (precision + recall))
            auc_roc_score = calculate_auc_roc(links, golden_links)
            label = DE_type + f' (score = )'
            ax_pr.plot(recall, precision, label=label, color=colors[i], alpha=1, linewidth=2,
                linestyle=linestyles[i])
            label = DE_type + f' (score = {str(round(auc_roc_score, 2))})'
            ax_roc.plot(fpr, tpr, label=label, color=colors[i], alpha=1, linewidth=2,
                       linestyle=linestyles[i])
            label = DE_type + f' (score = {str(round(np.mean(recall[:-1]), 2))})'
            ax.plot(thresholds, recall[:-1], label=label, color=colors[i], alpha=1, linewidth=2,
                        linestyle=linestyles[i])

        ax_pr.set_title(method)
        ax_pr.legend(frameon=False)
        ax_pr.set_ylabel('Precision')
        ax_pr.set_xlabel('Recall')

        ax_roc.set_title(method)
        ax_roc.legend(frameon=False)
        ax_roc.set_ylabel('True positive rate')
        ax_roc.set_xlabel('False positive rate')

        ax.set_title(method)
        ax.legend(frameon=False)
        ax.set_ylabel('Recall')
        ax.set_xlabel('Thresholds')

    # fig_pr.savefig(os.path.join(MODELSELECTION_DIR, f'precision_recall.png'), dpi=300, transparent=True)
    fig_pr.savefig(os.path.join(MODELSELECTION_DIR, f'precision_recall.pdf'))
    fig_roc.savefig(os.path.join(MODELSELECTION_DIR, f'roc.pdf'))
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'recall.pdf'))




