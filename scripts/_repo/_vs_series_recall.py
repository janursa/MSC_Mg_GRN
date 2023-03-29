"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing
import random


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, F_DE_protnames, CALIBRATION_DIR
from scripts.utils.links import plot_match_counts_series, normalize_links, format_links_string
from scripts.utils import calibration, serif_font, make_title_pretty

from geneRNI.evaluation import precision_recall_curve,calculate_auc_roc, calculate_PR


if __name__ == '__main__':
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)
    n_repeat = 10

    arbitrary_dir = os.path.join(GRN_DIR, 'arbitrary')
    if not os.path.isdir(arbitrary_dir):
        os.makedirs(arbitrary_dir)

    DE_types = F_DE_data().keys()
    study='all-in'
    # methods = ['RF', 'ridge', 'portia', 'arbitrary']
    methods = ['portia', 'arbitrary']

    # ncols = 2
    # nrows = 3
    # fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))

    links_string_dict = {}
    for DE_type in DE_types:
        links_string = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        links_string_dict[DE_type] = links_string

    # serif_font()
    # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey', 'lightcoral']
    # linestyles = ['-', '--', '-.', ':', '-', '--']

    for idx, DE_type in enumerate(DE_types):
        # i = int(idx/ncols)
        # j = idx%ncols
        # ax = axes[i][j]
        # - retreive or create the links
        links_stack = []

        for method in methods:
            if method in ['portia', 'ridge']:
                # - single weight
                links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
                links_stack.append(links)
            elif method == 'RF':
                links_pool = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
                links_stack.append(links_pool)
            elif method == 'arbitrary':
                links, links_pool = create_random_links(links_stack, n=n_repeat)
                links.to_csv(os.path.join(arbitrary_dir, f'links_{DE_type}_{study}.csv'), index=False)
                links_pool.to_pickle(os.path.join(arbitrary_dir, f'links_pool_{DE_type}_{study}.csv'))
                links_stack.append(links_pool)
        assert (len(links_stack)==2)
        #-- compare to golden links
        golden_links = format_links_string(links_string_dict[DE_type], F_DE_protnames()[DE_type])
        protnames = F_DE_protnames()[DE_type]
        #- remove mg from links
        def remove_mg(links):
            return links.loc[(links['Regulator']!='mg') & (links['Target']!='mg'),:]
        # - get the links, filter mask based on scores, and calculate recall
        for method_i, method in enumerate(methods):

            links = links_stack[method_i]
            links = remove_mg(links)
            # if method in ['ridge', 'RF']: #skipped
            #     # best_scores, best_params = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
            #     # protnames_left = np.asarray(protnames)[np.asarray(best_scores > 0)]
            #     # filter_mask = links['Target'].isin(protnames_left)
            #     filter_mask = None
            # else:
            #     filter_mask = None

            _, recall, _ = precision_recall_curve(links, golden_links)
            # calculate_auc_roc(links, golden_links)

            # label = method + f' (score = {str(round(np.mean(recall[:-1]), 2))})'
            # def normalize(aa):
            #     return (aa
            # ax.plot(, recall, label=label, color=colors[i], alpha=1, linewidth=2,
            #             linestyle=linestyles[i])
            # print(f'{DE_type} -- {method}---{round(np.mean(recall[:-1]), 2)}')
            print(f'auc roc -> {DE_type} -- {method}---{calculate_auc_roc(links, golden_links)}')
            print(f'PR -> {DE_type} -- {method}---{calculate_PR(links, golden_links)}')

        # ax.set_title(DE_type)
        # ax.legend(frameon=False)
        # ax.set_ylabel('Recall')
        # ax.set_xlabel('Thresholds')

    # fig_pr.savefig(os.path.join(MODELSELECTION_DIR, f'precision_recall.png'), dpi=300, transparent=True)
    # fig.savefig(os.path.join(MODELSELECTION_DIR, f'recall.pdf'))




