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

import tqdm

from scripts.imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, F_DE_protnames, CALIBRATION_DIR
from scripts.utils.links import choose_top_quantile, normalize_links
from scripts.utils import calibration, serif_font, make_title_pretty

from geneRNI.evaluation import precision_recall_curve,calculate_auc_roc, calculate_PR
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


def plot_dist(ax, data_stack, x_labels, sig_signs):
    matplotlib.rcParams.update({'font.size': 12})
    bplot = ax.violinplot(data_stack, showmeans=True, showextrema=False)
    ax.set_ylabel('Early precision ratio')
    ax.set_xticks(list(range(1,len(x_labels)+1)))
    ax.set_xticklabels(x_labels,rotation=0)
    ax.set_ymargin(.25)
    #- face colors
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan','grey']
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
def plot_line(ax, x_data, data_stack, line_names):
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey', 'lightcoral']
    linestyles = ['-', '--', '-.', ':', '-', '--']
    for i, line_name in enumerate(line_names):
        ax.plot(x_data,
                       data_stack[i],
                       label=line_names[i],
                       color=colors[i],
                       alpha=1,
                       linewidth=2,
                       linestyle=linestyles[i],
                       marker='o')

    ax_series.legend(frameon=False)
    ax_series.set_ylabel('Early precision ratio')
    ax_series.set_xlabel('Top quantile')
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
def create_random_links(links_assembly):
    #- TODO: create matrix (using protname)
    links_assembly = flatten(links_assembly)
    links_assembly = [normalize_links(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j] #flatten
    sample_links = links_assembly[0]
    random_links = pd.DataFrame({key:sample_links[key] for key in ['Regulator', 'Target']})
    random_links['Weight'] = random.sample(weights, len(random_links))
    return random_links

def calculate_EP(links, string_df, top_quantiles) -> int:
    '''  calculate early precision
        Compare extracted links by GRN to those suggested by model_selection.
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
        ns.append(n/len(golden_links))
    return np.asarray(ns)

if __name__ == '__main__':
    n_repeat = 40
    top_quantiles = np.linspace(.75, .9, 10)
    DE_types = F_DE_data().keys()
    study='ctr'
    methods = ['RF', 'ridge', 'portia', 'arbitrary']
    methods_preferred_names = ['RF', 'Ridge', 'Portia',  'Arbitrary'] # for line plot
    methods_preferred_order = ['Arbitrary', 'RF', 'Ridge', 'Portia'] # for dist plot

    #- plot specs
    ncols = 2
    nrows = 2
    fig_series, axes_series = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))
    fig_dist, axes_dist = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))

    #- get golden links for each DE_type
    links_string_dict = {}
    for DE_type in DE_types:
        links_string_dict[DE_type] = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
    #- create directories
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)

    #- main
    for idx, DE_type in enumerate(DE_types):
        golden_links = links_string_dict[DE_type]
        protnames = F_DE_protnames()[DE_type]
        # - retreive the links and create the random links
        links_stack = {}
        protnames_fs = [] #left protein names after filter based on score
        for method in methods:
            if method in ['portia', 'ridge', 'RF']:
                # - single weight
                links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
                protnames_f = None
            elif method == 'arbitrary':
                links: List[pd.DataFrame] = [] # it contains a list of links instead of a single links
                for i in range(n_repeat):
                    links.append(create_random_links(links_stack.values()))
                protnames_f = None
            links_stack[method] = links
            protnames_fs.append(protnames_f)

        #- calculate EP scores
        scores_series_stack: List[float] = [] #contains scores for line plot
        scores_dist_stack: List[float] = [] #contains scores for dist plot. The average score of each line
        for method_i, method in enumerate(methods):
            links = links_stack[method]
            if method in ['arbitrary']:
                # - multiple scores due to randomness
                scores_list: List[List[int]] = []
                for i in tqdm.tqdm(range(len(links)), desc=f'Progress {method} for {DE_type}'):
                    scores = calculate_EP(links[i], golden_links, top_quantiles)
                    scores_list.append(scores)
                scores_series = np.mean(scores_list, axis=0)
                scores_dist = np.mean(scores_list, axis=1)
            # elif method  in ['ridge_f','RF_f']:
            #     # - single score as there is no randomness
            #     protnames_f = protnames_fs[method_i]
            #     golden_links_f = golden_links.loc[golden_links['Target'].isin(protnames_f), :]
            #     scores_series = calculate_EP(links, golden_links_f, top_quantiles)
            #     scores_dist = np.repeat(np.mean(scores_series), n_repeat)
            #     if method == 'ridge_f':
            #         print(np.mean(scores_series))
            else:
                # - single score as there is no randomness
                scores_series = calculate_EP(links, golden_links, top_quantiles)
                # scores_dist = np.repeat(np.mean(scores_series), n_repeat)
                scores_dist = np.mean(scores_series)
            scores_series_stack.append(scores_series)
            scores_dist_stack.append(scores_dist)
        # try: # assert length equality
        #     np.vstack(scores_dist_stack)
        # except:
        #     raise ValueError('Probably length of scores do not match')

        # scores_series_stack = scores_series_stack/np.mean(scores_series_stack[-1]) #normalize with random
        # scores_dist_stack = scores_dist_stack/np.mean(scores_dist_stack[-1])

        # #- collect case spesific scores
        # for method_i, method in enumerate(methods):
        #     if method == 'arbitrary':
        #         continue
        #     method_scores[method].append(scores_dist_stack[method_i])
        # def sum_score(case):
        #     for key in case.keys():
        #         if key in DE_type:
        #             case[key].append(np.sum(scores_dist_stack[:-1]))
        # sum_score(imput_scores)
        # sum_score(phase_scores)
        # sum_score(pvalue_scores)
        #- line plot
        i = int(idx / ncols)
        j = idx % ncols
        ax_series = axes_series[i][j]
        plot_line(ax = ax_series, x_data=top_quantiles, data_stack=scores_series_stack, line_names=methods_preferred_names)
        title = DE_type+f'(prots: {len(protnames)}, slinks:{len(golden_links)})'
        ax_series.set_title(title)
        #- plot dist
        ax_dist = axes_dist[i][j]
        scores_dist_stack = list(scores_dist_stack)
        scores_dist_stack.insert(0, scores_dist_stack[-1]) #change to the preferred order
        del scores_dist_stack[-1]
        # - adjust
        scores_dist_stack_a = [scores_dist_stack[0]]
        for item in scores_dist_stack[1:]:  # except random
            scores_dist_stack_a.append(np.repeat(item, len(scores_dist_stack[0])))
        sig_signs = determine_sig_signes(scores_dist_stack_a)
        plot_dist(ax=ax_dist, data_stack=scores_dist_stack_a, x_labels=methods_preferred_order, sig_signs=sig_signs)
        ax_dist.set_title(DE_type)


    fig_series.savefig(os.path.join(MODELSELECTION_DIR, f'ep_vs_quantile.png'), dpi=300, transparent=True)
    fig_series.savefig(os.path.join(MODELSELECTION_DIR, f'ep_vs_quantile.pdf'))

    fig_dist.savefig(os.path.join(MODELSELECTION_DIR, f'ep_dist.png'), dpi=300, transparent=True)
    fig_dist.savefig(os.path.join(MODELSELECTION_DIR, f'ep_dist.pdf'))

    #-
    # def analyse(case):
    #     for key, item in case.items():
    #         print(f'key: {key}: {np.mean(item)}')
    # analyse(imput_scores)
    # analyse(phase_scores)
    # analyse(pvalue_scores)
    # analyse(method_scores)



