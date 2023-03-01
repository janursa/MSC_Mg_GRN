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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
sys.path.append(os.path.join(os.path.dirname(SCRIPT_DIR), '..'))
sys.path.append(os.path.join(os.path.dirname(SCRIPT_DIR), '..', '..'))

from typing import Dict, List
from scripts.utils import make_title_pretty

import tqdm

from scripts.imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, F_DE_protiens, CALIBRATION_DIR
from scripts.utils.links import choose_top_quantile, normalize_links
from scripts.utils import calibration, serif_font

def determine_sig_signes(data_stack):
    """ to test sig distribtion from noise links
    Conducts t test to determine whether datas[1:] are significantly different than datas[0], which is ctr
    Datas: Tuple(DataFrame), e.g. [ctr, RF, Ridge, Portia]
    """
    ctr = data_stack[0] #random
    #- determine p values: compared to ctr
    signs = ['']
    for data in data_stack[1:]:
        s, p = scipy.stats.ttest_ind(data, ctr)
        if p<0.05 and (np.mean(data) > np.mean(ctr)):
            sign = r'$*$'
        else:
            sign = ''
        signs.append(sign)
    return signs
def plot_dist(ax, idx, data_stack, x_labels, sig_signs, title):
    serif_font()
    # matplotlib.rcParams.update({'font.size': 12})
    bplot = ax.violinplot(data_stack.T, showmeans=True, showextrema=False)
    if idx %2 == 0:
        ax.set_ylabel('Early precision ratio')
        ax.set_yticks([0,1,2])
        ax.set_yticklabels([0, 1, 2])
    else:
        ax.set_yticks([0, 1, 2])
        ax.set_yticks([])
        ax.set_ylabel('')
    ax.set_xticks(list(range(1,len(x_labels)+1)))
    ax.set_xticklabels(x_labels,rotation=0)
    ax.set_ymargin(.25)
    #- face colors
    colors = ['lightgreen', 'lightblue', 'lightgreen', 'grey']
    for patch, color in zip(bplot['bodies'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(1)
    #- plot sig
    xs = ax.get_xticks()
    ys = np.max(data_stack, axis=1)
    for i, sign in enumerate(sig_signs):
        if sign != '':
            ax.annotate(sign, xy=(xs[i],ys[i]),
                ha='center',
                va='bottom',
                )
    ax.set_title(title)
def plot_line(ax, idx, x_data, data_stack, line_names, title):
    serif_font()
    # matplotlib.rcParams.update({'font.size': 10})
    colors = ['grey', 'lightpink', 'lightblue', 'lightgreen', 'blue', 'orange', 'cyan']
    linestyles = np.repeat(['-', '--', '-.', ':'], 2)
    for i, line_name in enumerate(line_names):
        ax.plot(x_data,
                       data_stack[i],
                       label=line_names[i],
                       color=colors[i],
                       alpha=1,
                       linewidth=2,
                       linestyle=linestyles[i],
                       marker='o')

    # ax.legend(frameon=False)
    if idx %2 == 0:
        ax.set_ylabel('Early precision ratio')
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels([0, 1, 2])
    else:
        ax.set_yticks([0, 1, 2])
        ax.set_yticks([])
        ax.set_ylabel('')
    ax.set_xlabel('Top quantile')
    ax.set_title(title)
    # ax.set_ymargin(.1)
    # ax.set_xmargin(.1)
def flatten(lst):
    """
    Flattens a list that may contain nested lists.
    """
    flattened = []
    for item in lst:
        if isinstance(item, list):
            flattened.extend(flatten(item))
        else:
            flattened.append(item)
    return flattened
def create_random_links(links_stack, n):
    """
    Creates n number of random links using weights of given set of links
    """
    links_stack = flatten(links_stack)
    links_stack = [normalize_links(links) for links in links_stack]
    weights = [links['Weight'].values.tolist() for links in links_stack]
    weights = [i for j in weights for i in j] #flatten
    sample_links = links_stack[0]
    random_links = []
    for i in range(n):
        rep = sample_links.copy()
        rep.loc[:,'Weight'] =random.sample(weights, len(rep))
        random_links.append(rep)

    return random_links

def calculate_EP(links, string_df, top_quantiles) -> List[float]:
    '''  calculate normalized early precision scores for each top_quantile
        Compare extracted links by GRN to those suggested by model_selection.
        Normalized with length of golden link to make it comparible between cases with different lengths.
    '''
    ns = []
    for tpq in top_quantiles:
        links_top = choose_top_quantile(links, quantile=tpq)

        l_regs = links_top['Regulator'].to_numpy(str)
        l_targs = links_top['Target'].to_numpy(str)

        s_regs = string_df['Regulator'].to_numpy(str)
        s_targs = string_df['Target'].to_numpy(str)

        n = 0
        for reg, target in zip(s_regs, s_targs):
            if ((l_regs == reg) & (l_targs == target)).any():
                n+=1
        ns.append(n)
    return np.asarray(ns)/len(golden_links)

if __name__ == '__main__':
    n_repeat = 100
    top_quantiles = np.linspace(.75, .9, 10)
    DE_types = F_DE_data().keys()
    study='ctr'
    methods = ['RF', 'ridge', 'portia']

    #- plot specs
    ncols = 2
    nrows = 2
    fig_series, axes_series = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.8 * ncols, 2.4 * nrows))
    fig_dist, axes_dist = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2 * nrows))

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
        exclude_flags = [False for i in methods] # if True, the model will not be considered for further evaluation
        # - retreive the links
        links_stack = []
        for method in methods:
            if method in ['portia', 'ridge', 'RF']:
                links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
            links_stack.append(links)
        #-  create  random links
        links_random = create_random_links(links_stack, n_repeat)
        #- retreive the test scores
        test_scores = [] #if the mean best scores across both studies is greater than 0
        for method_i, method in enumerate(methods):
            if method in ['ridge', 'RF']:
                scores_studies = []
                for study in ['ctr', 'mg']:
                    best_scores, _ = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
                    mean_score = np.mean(best_scores)
                    if mean_score <= 0:
                        exclude_flags[method_i] = True
                    scores_studies.append(mean_score)
                test_score = np.mean(scores_studies)
            else:
                test_score = None
            test_scores.append(test_score)
        #- calculate early precision scores (normalized to the length of golden links)
        # random links
        golden_links = links_string_dict[DE_type]
        ep_scores_random: List[float] = []
        ep_scores_random_list: List[List[float]] = [] # for line plot
        for i in tqdm.tqdm(range(n_repeat), desc=f'Progress for {DE_type}'):
            scores = calculate_EP(links_random[i], golden_links, top_quantiles)
            ep_scores_random_list.append(scores)
            ep_scores_random.append(np.mean(scores))
        # other methods
        ep_scores: List[float] = []
        ep_scores_list: List[List[float]] = [] # for line plot
        for method_i,_ in enumerate(methods):
            links = links_stack[method_i]
            scores = calculate_EP(links, golden_links, top_quantiles)
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
        protnames = F_DE_protiens()[DE_type]
        methods_preferred_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']
        # line plot
        i = int(idx / ncols)
        j = idx % ncols
        ax_series = axes_series[i][j]
        epr_scores_list_random_included = np.vstack([[1 for _ in top_quantiles], epr_scores_list])
        plot_line(ax = ax_series, idx=idx, x_data=top_quantiles, data_stack=epr_scores_list_random_included, line_names=methods_preferred_names, title=make_title_pretty(DE_type))
        if idx == 1:
            ax_series.legend(loc='upper center', bbox_to_anchor=(1.35, 1), fancybox=False, frameon=False)

        # title = DE_type+f'(prots: {len(protnames)}, slinks:{len(golden_links)})'
        # dist plot
        ax_dist = axes_dist[i][j]
        epr_scores_random = ep_scores_random / np.mean(ep_scores_random)
        sig_signs = [r'$*$' if (p < 0.05) & (score>1) else '' for p,score in zip(pvalues,epr_scores)]
        sig_signs.insert(0, '')  # for random
        scores_dist = np.vstack([epr_scores_random, [np.repeat(score, n_repeat) for score in epr_scores]])
        plot_dist(ax=ax_dist, idx=idx, data_stack=scores_dist, x_labels=methods_preferred_names, sig_signs=sig_signs, title=make_title_pretty(DE_type))

    fig_series.savefig(os.path.join(MODELSELECTION_DIR, f'ep_vs_quantile.png'), dpi=300, transparent=True)
    fig_series.savefig(os.path.join(MODELSELECTION_DIR, f'ep_vs_quantile.pdf'))

    fig_dist.savefig(os.path.join(MODELSELECTION_DIR, f'ep_dist.png'), dpi=300, transparent=True)
    fig_dist.savefig(os.path.join(MODELSELECTION_DIR, f'ep_dist.pdf'))

    #-plot the main graph for the selected models
    model_names_with_scores = [f'{name} ({round(np.mean(score),2)})' for name, score in zip(shortlisted_models_names,shortlisted_models_scores)]
    model_names_randome_included = [f'Arbitrary (1)']+ model_names_with_scores
    scores_random_included = np.vstack([[1 for _ in top_quantiles], shortlisted_models_scores])

    fig_main, ax_main = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.5))
    plot_line(ax=ax_main, idx=0, x_data=top_quantiles, data_stack=scores_random_included,
              line_names=model_names_randome_included, title='')

    ax_main.legend(loc='upper center', bbox_to_anchor=(1.48, 1), fancybox=False, frameon=False)
    fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'main.png'), bbox_inches="tight", dpi=300, transparent=True)
    fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'main.pdf'), bbox_inches="tight")

    # np.savetxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), shortlisted_models_names, delimiter=",", fmt="%s")






