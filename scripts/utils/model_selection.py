"""
    Sets of functions useful for model selection
"""

import sys
import os
import numpy as np
import random

import pandas as pd
import scipy
from typing import List, Dict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import serif_font, flatten
from utils.links import choose_top_quantile, normalize_links

def violinplot(ax, data_stack:np.ndarray, x_labels:List[str], sig_signs:List[str], percentages:List[str], title:str):
    # - x axis
    x_tick_locations = [i + 1 + (i // 8) for i in range(len(x_labels))] # every 8 groups closer to each other
    ax.set_xticks(x_tick_locations)
    ax.set_xticklabels(x_labels, rotation=30)
    ax.set_xmargin(.025)
    # - y axis
    ax.set_ymargin(.25)
    ax.set_ylabel('Mean EPR')
    # - violin plot
    serif_font()
    group_colors = ["limegreen", "orangered"]
    bplot = ax.violinplot(data_stack.T, positions=x_tick_locations, showmeans=False, showextrema=False)
    # - plot medians as scatter plot
    quartile1, medians, quartile3 = np.percentile(data_stack, [25, 50, 75], axis=1)
    group_repeated_colors = list(np.repeat(group_colors, 4, axis=0))
    meancolors = group_repeated_colors + group_repeated_colors
    ax.scatter(x_tick_locations, medians, marker='o', color=meancolors, s=30, zorder=3)

    #- face colors
    for patch_i, patch in enumerate(bplot['bodies']):
        patch.set_facecolor(meancolors[patch_i])
        patch.set_edgecolor('black')
        patch.set_alpha(1)
    #- plot sig
    xs = ax.get_xticks()
    ys = np.max(data_stack, axis=1)
    for i, value in enumerate(percentages):
        if value != '':
            ax.annotate(value, xy=(xs[i],ys[i]),
                ha='center',
                va='bottom',
                verticalalignment='baseline',
                textcoords='offset points',
                xytext=(0, 5)
                        )
    # - plot significannce
    for i, value in enumerate(sig_signs):
        if value != '':
            ax.annotate(value, xy=(xs[i], ys[i]),
                        ha='center',
                        va='bottom',
                        verticalalignment='baseline',
                        textcoords='offset points',
                        xytext=(0, -15)
                        )
    # - plot the legends for the groups
    handles = []
    labels = ['MinProb', 'KNN']
    for i, color in enumerate(group_colors):
        handles.append(ax.scatter([], [], marker='o', label=labels[i], color=color,
                                  s=30, alpha=1))
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, handles=handles, title='', fancybox=False,
                   frameon=False, prop={'size': 10}, title_fontproperties={'size': 9, 'weight': 'bold'})

    ax.set_title(title)
def lineplot(ax, idx, x_data, data_stack, line_names, title, yticks):
    """ Line plot for the epr scores versus top quantiles

    """
    serif_font()
    colors = ['grey', 'lightpink', 'lightblue', 'lightgreen', 'blue', 'orange', 'cyan', 'cyan', 'cyan','cyan']
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
    ax.set_ylabel('Early precision ratio')
    ax.set_ymargin(.2)
    ax.set_xlabel('Top quantile')
    ax.set_title(title)
    if yticks is not None:
        ax.set_yticks(yticks)

def create_random_links(links_stack: List[pd.DataFrame], n:int) -> List[pd.DataFrame]:
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
        rep.loc[:,'Weight'] = random.sample(weights, len(rep))
        random_links.append(rep)

    return random_links

def calculate_early_precision(links, golden_links, top_quantiles) -> List[float]:
    '''  calculate normalized EP scores for the top_quantiles
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
    return np.asarray(ns)/len(golden_links)*100
