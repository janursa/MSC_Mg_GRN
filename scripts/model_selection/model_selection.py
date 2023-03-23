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
import time
import matplotlib
import scipy
import tqdm
import argparse


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from typing import Dict, List
from utils import make_title_pretty


from imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, F_DE_proteins, CALIBRATION_DIR
from utils import calibration
from utils.model_selection import create_random_links, calculate_early_precision, violinplot, lineplot

def retreive_links(methods, DE_type, study):
    # - retreive the links
    links_stack = []
    for method in methods:
        links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
        links_stack.append(links)
    return links_stack
def retrieve_test_scores(methods, DE_type):
    test_scores = []  # if the mean best scores across both studies is greater than 0
    for method_i, method in enumerate(methods):
        if method in ['ridge', 'RF']:
            scores_studies = []
            for study in studies:
                best_scores, _ = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
                mean_score = np.mean(best_scores)
                # if mean_score <= 0:
                #     exclude_flags[method_i] = True
                scores_studies.append(mean_score)
            # test_score = np.mean(scores_studies)
        else:
            scores_studies = None
        test_scores.append(scores_studies)
    return test_scores
def filter_based_on_test_scores(test_scores):
    exclude_flags = [False for i in methods]  # if True, the model will not be considered for further evaluation
    for method_i, method_score in enumerate(test_scores):
        if method_score is None:
            continue
        for study_i, study in enumerate(studies):
            method_score
            study_score = method_score[study_i]
            if study_score <= 0:
                exclude_flags[method_i] = True
    return exclude_flags


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--methods', nargs='+', default=['RF', 'ridge', 'portia'])
    parser.add_argument('--n_repeat', nargs='+', default=100, help='Number of randomly generated noisy links')
    parser.add_argument('--top_quantiles', default=np.linspace(.75, .9, 10))
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    methods = args.methods
    n_repeat = args.n_repeat
    top_quantiles = args.top_quantiles

    DE_types = F_DE_data().keys()

    #- plot specs
    ncols = 2
    nrows = 2
    fig_series, axes_series = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(3 * ncols, 2.25 * nrows))
    fig_dist, axes_dist = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2.25 * nrows))

    #- get golden links for each DE_type
    links_string_dict = {}
    for DE_type in DE_types:
        links_string_dict[DE_type] = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
    #- create directories
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)

    #- main
    shortlisted_models_names = []
    shortlisted_models_scores = [] # remaining models
    for idx, DE_type in enumerate(DE_types):
        links_stack = retreive_links(methods, DE_type, study='ctr')
        #-  create  random links
        links_random = create_random_links(links_stack, n_repeat)
        #- retreive the test scores
        test_scores = retrieve_test_scores(methods, DE_type)
        #- filter based on test scores
        exclude_flags = filter_based_on_test_scores(test_scores)
        #- calculate early precision scores (normalized to the length of golden links)
        # random links
        golden_links = links_string_dict[DE_type]
        ep_scores_random: List[float] = []
        ep_scores_random_list: List[List[float]] = [] # for line plot
        for i in tqdm.tqdm(range(n_repeat), desc=f'Progress for {DE_type}'):
            scores = calculate_early_precision(links_random[i], golden_links, top_quantiles)
            ep_scores_random_list.append(scores)
            ep_scores_random.append(np.mean(scores))
        # other methods
        ep_scores: List[float] = []
        ep_scores_list: List[List[float]] = [] # for line plot
        for method_i,_ in enumerate(methods):
            links = links_stack[method_i]
            scores = calculate_early_precision(links, golden_links, top_quantiles)
            ep_scores_list.append(scores)
            ep_scores.append(np.mean(scores))
        # - exclude those with non-significant ep scores
        pvalues = []
        for method_i, _ in enumerate(methods):
            ep_score = ep_scores[method_i]
            s, p = scipy.stats.ttest_1samp(ep_scores_random,ep_score)
            pvalues.append(p)
            if p >= 0.05 or (ep_score <= np.mean(ep_scores_random)):
                exclude_flags[method_i] = True
        # - calculate early precision ratio
        epr_scores = ep_scores/np.mean(ep_scores_random)
        epr_scores_list = ep_scores_list/np.mean(ep_scores_random_list, axis=0) #this is divided by the mean value of line
        # epr_scores_random_list = ep_scores_random_list/np.mean(ep_scores_random_list)
        # - store selected models
        for method_i, method in enumerate(methods):
            if not exclude_flags[method_i]:
                tag = '_'.join([DE_type, method])
                # tag+=f' ({round(epr_scores[method_i],2)})'
                shortlisted_models_names.append(make_title_pretty(tag))
                shortlisted_models_scores.append(epr_scores_list[method_i])

        #- plots
        protnames = F_DE_proteins()[DE_type]
        methods_preferred_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']
        # line plot
        i = int(idx / ncols)
        j = idx % ncols
        ax_series = axes_series[i][j]
        epr_scores_list_random_included = np.vstack([[1 for _ in top_quantiles], epr_scores_list])
        lineplot(ax = ax_series, idx=idx, x_data=top_quantiles, data_stack=epr_scores_list_random_included, line_names=methods_preferred_names, title=make_title_pretty(DE_type), yticks=[0, 1, 2])
        if idx == 1:
            ax_series.legend(loc='upper center', bbox_to_anchor=(1.35, 1), fancybox=False, frameon=False)

        # title = DE_type+f'(prots: {len(protnames)}, slinks:{len(golden_links)})'
        # dist plot
        ax_dist = axes_dist[i][j]
        epr_scores_random = ep_scores_random / np.mean(ep_scores_random)
        sig_signs = [r'$*$' if (p < 0.05) & (score>1) else '' for p,score in zip(pvalues,epr_scores)]
        sig_signs.insert(0, '')  # for random
        scores_dist = np.vstack([epr_scores_random, [np.repeat(score, n_repeat) for score in epr_scores]])
        violinplot(ax=ax_dist, idx=idx, data_stack=scores_dist, x_labels=methods_preferred_names, sig_signs=sig_signs, title=make_title_pretty(DE_type))

    fig_series.savefig(os.path.join(MODELSELECTION_DIR, f'ep_vs_quantile.png'), dpi=300, transparent=True)
    fig_series.savefig(os.path.join(MODELSELECTION_DIR, f'ep_vs_quantile.pdf'))

    fig_dist.savefig(os.path.join(MODELSELECTION_DIR, f'ep_dist.png'), dpi=300, transparent=True)
    fig_dist.savefig(os.path.join(MODELSELECTION_DIR, f'ep_dist.pdf'))

    #-plot the main graph for the selected models
    model_names_with_scores = [f'{name} ({round(np.sum(score),2)})' for name, score in zip(shortlisted_models_names,shortlisted_models_scores)]
    model_names_randome_included = [f'Arbitrary ({len(top_quantiles)})']+ model_names_with_scores
    scores_random_included = np.vstack([[1 for _ in top_quantiles], shortlisted_models_scores])

    fig_main, ax_main = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.5))
    lineplot(ax=ax_main, idx=0, x_data=top_quantiles, data_stack=scores_random_included,
              line_names=model_names_randome_included, title='', yticks= None)

    ax_main.legend(loc='upper center', bbox_to_anchor=(1.48, 1), fancybox=False, frameon=False, title='Model (AOC)', title_fontproperties={'weight':'bold'} )
    fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'shortlisted.png'), bbox_inches="tight", dpi=300, transparent=True)
    fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'shortlisted.pdf'), bbox_inches="tight")

    # np.savetxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), shortlisted_models_names, delimiter=",", fmt="%s")






