"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import json
import argparse
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
from typing import Dict, List, Tuple, Callable, Optional, TypeAlias, Any



SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import MODELSELECTION_DIR, make_title_pretty
from md_aux import retreieve_scores
from utils import serif_font

def is_single_value(value):
    collection_types = (list, tuple, set, dict, np.ndarray)
    return not isinstance(value, collection_types)
def violinplot(data_stack: List[List[float]], percentages:List[str], sig_signs:List[str], GRN_methods:List[str]) -> Figure:
    """Plots violin for the distribution of epr scores"""
    fig, ax = plt.subplots(nrows=1, ncols=1, tight_layout=True, figsize=(6, 2.2))
    # - x axis
    x_tick_locations = [i + 1 + (i // 8) for i in range(len(GRN_methods))] # every 8 groups closer to each other
    ax.set_xticks(x_tick_locations)
    ax.set_xticklabels(GRN_methods, rotation=45, fontsize=9)
    ax.set_xmargin(.025)
    # - y axis
    ax.set_ymargin(.15)
    ax.set_ylabel('EPR (mean)')
    # - violin plot
    serif_font()
    group_colors = ["#1E90FF", "#FFA07A"]
    bplot = ax.violinplot(np.asarray(data_stack).T, positions=x_tick_locations, showmeans=False, showextrema=False)
    # - plot medians as scatter plot
    quartile1, medians, quartile3 = np.percentile(data_stack, [25, 50, 75], axis=1)
    group_repeated_colors = list(np.repeat(group_colors, 4, axis=0))
    meancolors = group_repeated_colors + group_repeated_colors
    ax.scatter(x_tick_locations, medians, marker='o', color=meancolors, s=45, zorder=3)
    #- face colors
    for patch_i, patch in enumerate(bplot['bodies']):
        patch.set_facecolor(meancolors[patch_i])
        patch.set_edgecolor('gray')
        patch.set_alpha(1)
    #- annotate percentile rank
    xs = ax.get_xticks()
    ys = np.max(data_stack, axis=1)
    for i, value in enumerate(percentages):
        if value != '':
            ax.annotate(value, xy=(xs[i],ys[i]),
                ha='center',
                va='bottom',
                verticalalignment='baseline',
                textcoords='offset points',
                fontsize = 9,
                xytext=(0, 5)
                        )
    # - annotate sig sign
    for i, value in enumerate(sig_signs):
        ax.annotate(value, xy=(xs[i], ys[i]),
                        ha='center',
                        va='bottom',
                        verticalalignment='baseline',
                        textcoords='offset points',
                        color='red',
                        xytext=(0, -15)
                        )
    # # - plot the legends for the groups
    # handles = []
    # labels = ['MinProb', 'KNN']
    # for i, color in enumerate(group_colors):
    #     handles.append(ax.scatter([], [], marker='o', label=labels[i], color=color,
    #                               s=30, alpha=1))
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, handles=handles, title='', fancybox=False,
    #                frameon=False, prop={'size': 10}, title_fontproperties={'size': 9, 'weight': 'bold'})

    return fig
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_random_links', default=1000, type=int, help='Number of randomly generated noisy links')
    args, remaining_args = parser.parse_known_args()

    n_repeat = args.n_random_links

    #- retreieve scores
    scores = retreieve_scores()
    # - plot violin
    methods_names = [make_title_pretty(score.split('_')[2]) for score in scores.keys()] # only method name
    scores_list = [item['epr'] for item in scores.values()]
    sig_flags = [item['sig_flag'] for item in scores.values()]
    percentile_ranks = [item['percentile_rank'] for item in scores.values()]

    sig_signs = [r'$*$' if flag else '' for flag in sig_flags]
    percentile_ranks = [f'{item}%' if item else '' for item in percentile_ranks]

    scores_dist = [np.repeat(score, n_repeat) if (is_single_value(score)) else score for score in scores_list]

    fig = violinplot(scores_dist, percentile_ranks, sig_signs, methods_names)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot.pdf'))

