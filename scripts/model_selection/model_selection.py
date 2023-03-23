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

def retrieve_test_scores(DE_type, method, studies):
    if method == 'portia':
        scores_studies = None
    else:
        scores_studies = []
        for study in studies:
            best_scores, _ = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
            mean_score = np.mean(best_scores)
            scores_studies.append(mean_score)
    return scores_studies
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

def calculate_early_precision_ratio(links_methods, links_random, DE_type, top_quantiles):
    golden_links = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
    roc_scores_random: List[float] = []
    ep_scores_random_series: List[List[float]] = []  # 100*10
    # - epr scores for random
    for i in tqdm.tqdm(range(len(links_random)), desc=f'Progress for {DE_type}'):
        scores = calculate_early_precision(links_random[i], golden_links, top_quantiles) #10 scores one for each top quantile
        ep_scores_random_series.append(scores)
        roc_scores_random.append(np.sum(scores))
    # - epr scores for the rest of methods
    roc_scores: List[float] = []
    ep_scores_series: List[List[float]] = []  # 100*10
    for method_i,method in enumerate(methods):
        scores = calculate_early_precision(links_methods[method_i], golden_links, top_quantiles)
        ep_scores_series.append(scores)
        roc_scores.append(np.sum(scores))

    roc_scores_n = roc_scores / np.mean(roc_scores_random)
    roc_scores_random_n = roc_scores_random/np.mean(roc_scores_random)
    epr_scores_series = ep_scores_series / np.mean(ep_scores_random_series,
                                                             axis=0)  # this is divided by the mean value of line

    return roc_scores_n, roc_scores_random_n, epr_scores_series

def retreive_links(DE_types, methods):
    method_links_stack = {}  # only for ctr -> we need it for epr calculation
    random_links_stack = {}
    for DE_type in DE_types:
        links_methods = []
        for method in methods:
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_ctr.csv'), index_col=False)
            links_methods.append(links)
        random_links = create_random_links(links_methods, n_repeat)

        method_links_stack[DE_type] = links_methods
        random_links_stack[DE_type] = random_links
    return method_links_stack, random_links_stack
def extract_tag_based_data_from_modeldata(data, tag):
    selected = []
    for key, value in data.items():
        if tag in key:
            selected.append(value)
    return selected
def lineplot_all_models(epr_scores_series, top_quantiles):
    ncols = 2
    nrows = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True,
                                           figsize=(3 * ncols, 2.25 * nrows))

    for idx, DE_type in enumerate(DE_types):
        i = int(idx / ncols)
        j = idx % ncols
        ax_series = axes[i][j]
        eprs = extract_tag_based_data_from_modeldata(epr_scores_series, DE_type)
        epr_scores_all = np.vstack([[1 for _ in top_quantiles], eprs])
        lineplot(ax=ax_series, idx=idx, x_data=top_quantiles, data_stack=epr_scores_all,
                 line_names=methods_preferred_names, title=make_title_pretty(DE_type), yticks=[0, 1, 2])
        if idx == 1:
            ax_series.legend(loc='upper center', bbox_to_anchor=(1.35, 1), fancybox=False, frameon=False)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_all_models.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_all_models.pdf'))
def determine_p_value(epr_scores_methods, random_values):
    pvalues = []
    for method_i, _ in enumerate(methods):
        epr_score = epr_scores_methods[method_i]
        s, p = scipy.stats.ttest_1samp(random_values, epr_score)
        pvalues.append(p)
    return pvalues
def violinplot_all_models(roc_scores, roc_scores_random):

            # if p >= 0.05 or (ep_score <= np.mean(ep_scores_random)):
            #     exclude_flags[method_i] = True

    ncols = 2
    nrows = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2.25 * nrows))

    for idx, DE_type in enumerate(DE_types):
        random_values = roc_scores_random[DE_type]
        roc_scores_methods = extract_tag_based_data_from_modeldata(roc_scores, DE_type)
        pvalues = determine_p_value(roc_scores_methods, random_values)
        i = int(idx / ncols)
        j = idx % ncols
        ax_dist = axes[i][j]
        sig_signs = [r'$*$' if (p < 0.05) & (score>1) else '' for p, score in zip(pvalues, roc_scores_methods)]
        sig_signs.insert(0, '')  # for random
        scores_dist = np.vstack([random_values, [np.repeat(score, n_repeat) for score in roc_scores_methods]])
        violinplot(ax=ax_dist, idx=idx, data_stack=scores_dist, x_labels=methods_preferred_names, sig_signs=sig_signs, title=make_title_pretty(DE_type))

    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot_all_models.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot_all_models.pdf'))

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

    methods_preferred_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']

    DE_types = F_DE_data().keys()

    # - create directories
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)

    #- retreive links and calculate random links
    links_methods_stack, links_random_stack = retreive_links(DE_types, methods)
    #- calculate test score and epr for each model
    test_scores = {} # all models
    roc_scores = {} # all models
    epr_scores_series = {} # for line plot of epr vs top quantiles. DE_type: scores_series_methods
    roc_scores_random = {}
    for idx, DE_type in enumerate(DE_types):
        roc_scores_methods, roc_scores_random_methods, epr_scores_series_methods = calculate_early_precision_ratio(links_methods=links_methods_stack[DE_type],links_random=links_random_stack[DE_type], DE_type=DE_type, top_quantiles=top_quantiles)
        roc_scores_random[DE_type] = roc_scores_random_methods
        for method_i, method in enumerate(methods):
            model_name = '_'.join([DE_type, method])

            test_score = retrieve_test_scores(DE_type, method, studies) # test score is calculated for both ctr and sample

            test_scores[model_name] = test_score
            roc_scores[model_name] = roc_scores_methods[method_i]
            epr_scores_series[model_name] = epr_scores_series_methods[method_i]
            # epr_scores_series[model_name] = epr_scores_series
    lineplot_all_models(epr_scores_series, top_quantiles)
    violinplot_all_models(roc_scores, roc_scores_random)
    aa
    #- filter models with test_score less than the baseline model && those with epr_score less than random links
    # shortlisted_models = []
    # for model_name, scores_model in scores.item():
    #     scores_model
    #         links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
    # #- model selection
    # shortlisted_models_names = []
    # shortlisted_models_scores = []  # remaining models
    #
    #     #- filter based on test scores
    #     exclude_flags = filter_based_on_test_scores(test_scores)
    #
    #     # - exclude those with non-significant ep scores
    #     pvalues = []
    #     for method_i, _ in enumerate(methods):
    #         ep_score = ep_scores[method_i]
    #         s, p = scipy.stats.ttest_1samp(ep_scores_random,ep_score)
    #         pvalues.append(p)
    #         if p >= 0.05 or (ep_score <= np.mean(ep_scores_random)):
    #             exclude_flags[method_i] = True
    #     # epr_scores_random_list = ep_scores_random_list/np.mean(ep_scores_random_list)
    #     # - store selected models
    #     for method_i, method in enumerate(methods):
    #         if not exclude_flags[method_i]:
    #             tag = '_'.join([DE_type, method])
    #             # tag+=f' ({round(epr_scores[method_i],2)})'
    #             shortlisted_models_names.append(make_title_pretty(tag))
    #             shortlisted_models_scores.append(epr_scores_list[method_i])
    #
    #     #- plots
    #
    # #-plot the main graph for the selected models
    # model_names_with_scores = [f'{name} ({round(np.sum(score),2)})' for name, score in zip(shortlisted_models_names,shortlisted_models_scores)]
    # model_names_randome_included = [f'Arbitrary ({len(top_quantiles)})']+ model_names_with_scores
    # scores_random_included = np.vstack([[1 for _ in top_quantiles], shortlisted_models_scores])
    #
    # fig_main, ax_main = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.5))
    # lineplot(ax=ax_main, idx=0, x_data=top_quantiles, data_stack=scores_random_included,
    #           line_names=model_names_randome_included, title='', yticks= None)
    #
    # ax_main.legend(loc='upper center', bbox_to_anchor=(1.48, 1), fancybox=False, frameon=False, title='Model (AOC)', title_fontproperties={'weight':'bold'} )
    # fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'shortlisted.png'), bbox_inches="tight", dpi=300, transparent=True)
    # fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'shortlisted.pdf'), bbox_inches="tight")

    # np.savetxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), shortlisted_models_names, delimiter=",", fmt="%s")
    #
    #
    #
    #


