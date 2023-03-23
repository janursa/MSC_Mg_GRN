"""
    Sets of functions useful for hyperparameter tuning
"""
import sys
import os
import matplotlib.pylab as plt
import numpy as np
import random
from typing import List, Dict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import serif_font
from utils.links import choose_top_quantile, normalize_links, flatten
def violinplot(ax, idx, data_stack, x_labels, sig_signs, title):
    serif_font()
    # matplotlib.rcParams.update({'font.size': 12})
    bplot = ax.violinplot(data_stack.T, showmeans=True, showextrema=False)
    if idx %2 == 0:
        ax.set_ylabel('AOC-EPR (standardize)')
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
def lineplot(ax, idx, x_data, data_stack, line_names, title, yticks):
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
    # if idx %2 == 0:
    #     ax.set_ylabel('EPR')
    #     ax.set_yticks(yticks)
    #     ax.set_yticklabels(yticks)
    # else:
    #     ax.set_yticks(yticks)
    #     ax.set_yticks([])
    #     ax.set_ylabel('')
    ax.set_ylabel('Early precision ratio')
    ax.set_ymargin(.2)
    # ax.set_yticks(yticks)
    # ax.set_yticklabels(yticks)

    ax.set_xlabel('Top quantile')
    # if idx in [0,1]:
    #     ax.set_xlabel('')
    #     # ax.set_xticks([])
    #     # ax.set_xticklabels([])


    ax.set_title(title)
    # ax.set_ymargin(.1)
    # ax.set_xmargin(.1)

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

def calculate_early_precision(links, golden_links, top_quantiles) -> List[float]:
    '''  calculate normalized early precision scores for the top_quantiles
        Compare extracted links by GRN to those suggested by model_selection.
        Normalized with length of golden link to make it comparible between cases with different lengths.
    '''
    ns = []
    for tpq in top_quantiles:
        links_top = choose_top_quantile(links, quantile=tpq)

        l_regs = links_top['Regulator'].to_numpy(str)
        l_targs = links_top['Target'].to_numpy(str)

        s_regs = golden_links['Regulator'].to_numpy(str)
        s_targs = golden_links['Target'].to_numpy(str)

        n = 0
        for reg, target in zip(s_regs, s_targs):
            if ((l_regs == reg) & (l_targs == target)).any():
                n+=1
        ns.append(n)
    return np.asarray(ns)/len(golden_links)