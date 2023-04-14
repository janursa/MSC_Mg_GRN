"""
    Auxillary functions for model selection
"""
import sys
import os
import numpy as np
import json
from pathlib import Path
from typing import List, TypeAlias, Dict, Any, Tuple, Optional
import random
import pandas as pd
from typing import List, Dict
import scipy

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import MODELSELECTION_DIR, RANDOM_MODELS_DIR
from utils import serif_font, flatten
from utils.links import choose_top_quantile, normalize_links

score_type: TypeAlias = Dict[str, Any]  # template to store score. score_name, e.g. ep: score value


def save_scores(scores: Dict[str, score_type]) -> None:
    with open(f'{MODELSELECTION_DIR}/scores.json', 'w') as f:
        json.dump(scores, f)


def retreieve_scores() -> Dict[str, score_type]:
    with open(f'{MODELSELECTION_DIR}/scores.json', 'r') as f:
        # load the JSON data from the file
        scores = json.load(f)
    return scores


def get_baseline_scores(data_type: str) -> Tuple[List[float], List[List[float]]]:
    """Retreive baseline (random) scores of ep-mean"""
    ep_scores = np.loadtxt(Path(RANDOM_MODELS_DIR) / f'ep_scores_{data_type}.csv', delimiter=',')
    ep_scores_series = np.loadtxt(Path(RANDOM_MODELS_DIR) / f'ep_scores_series_{data_type}.csv', delimiter=',')

    return list(ep_scores), list(ep_scores_series)


def save_baseline_scores(ep_scores: List[float], ep_scores_series: List[List[float]], DE_type: str) -> None:
    np.savetxt(f'{Path(RANDOM_MODELS_DIR)}/ep_scores_{DE_type}.csv', np.asarray(ep_scores),
               delimiter=',')
    np.savetxt(f'{Path(RANDOM_MODELS_DIR)}/ep_scores_series_{DE_type}.csv', np.asarray(ep_scores_series),
               delimiter=',')


def lineplot(ax, x_data, data_stack, line_names, title, yticks):
    """ Line plot for the epr scores versus top quantiles

    """
    serif_font()
    colors = ['grey', 'lightpink', 'lightblue', 'lightgreen', 'blue', 'orange', 'cyan', 'cyan', 'cyan', 'cyan']
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


def create_random_links(links_stack: List[pd.DataFrame], n: int) -> List[pd.DataFrame]:
    """
    Creates n number of random links using weights of given set of links
    """
    links_stack = flatten(links_stack)
    links_stack = [normalize_links(links) for links in links_stack]
    weights = [links['Weight'].values.tolist() for links in links_stack]
    weights = [i for j in weights for i in j]  # flatten
    sample_links = links_stack[0]
    random_links = []
    for i in range(n):
        rep = sample_links.copy()
        rep.loc[:, 'Weight'] = random.sample(weights, len(rep))
        random_links.append(rep)

    return random_links


def is_single_value(value: Any) -> bool:
    """check if the given value is a single value"""
    collection_types = (list, tuple, set, dict, np.ndarray)
    return not isinstance(value, collection_types)


def determine_sig_flag(random_scores, score):
    s, p = scipy.stats.ttest_1samp(random_scores, score)
    if (p < 0.05) & (score > np.mean(random_scores)):
        sig_flag = True
    else:
        sig_flag = False
    return sig_flag


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
                n += 1
        ns.append(n)
    return list(np.asarray(ns) / len(golden_links) * 100)


class ViolinPlot:
    group_colors = ["#1E90FF", "#FFA07A"]

    @staticmethod
    def modify_x_axis(ax, ticks, ticknames, x_margin):
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticknames, rotation=45, fontsize=9)
        if x_margin is None:
            x_margin = 0.025
        ax.set_xmargin(x_margin)

    @staticmethod
    def modify_y_axis(ax, label):
        ax.set_ymargin(.15)
        ax.set_ylabel(label)

    @staticmethod
    def annotate_percentile_rank(ax, xs: List[float], ys: List[float], percentages: List[float]) -> None:
        """Write percentile ranks on top of each scatter point"""
        percentages = [f'{item}%' if item else '' for item in percentages]

        for i, value in enumerate(percentages):
            if value != '':
                ax.annotate(value, xy=(xs[i], ys[i]),
                            ha='center',
                            va='bottom',
                            verticalalignment='baseline',
                            textcoords='offset points',
                            fontsize=9,
                            xytext=(0, 5)
                            )

    @staticmethod
    def annotate_sig_signs(ax, xs: List[float], ys: List[float], sig_flags: List[bool], colors: List[str]) -> None:
        """Writes significant signs (*) on buttom of each scatter point"""
        sig_signs = [r'$*$' if flag else '' for flag in sig_flags]
        for i, value in enumerate(sig_signs):
            ax.annotate(value, xy=(xs[i], ys[i]),
                        ha='center',
                        va='bottom',
                        verticalalignment='baseline',
                        textcoords='offset points',
                        size=14,
                        color=colors[i],
                        xytext=(0, -20)
                        )

    @classmethod
    def place_legends_for_groups(cls, ax, legends) -> None:
        if not legends[0]:
            return
        handles = []
        if legends[1] is None:
            legends[1] = ['MinProb', 'KNN']
        for i, color in enumerate(cls.group_colors):
            handles.append(ax.scatter([], [], marker='o', label=legends[1][i], color=color,
                                      s=30, alpha=1))
        ax.legend(loc='upper left', bbox_to_anchor=(.5, .5), ncol=1, handles=handles, title='', fancybox=False,
                  frameon=False, prop={'size': 10}, title_fontproperties={'size': 9, 'weight': 'bold'})

    @staticmethod
    def plot_median_as_scatter_points(ax, xs: List[float], ys: List[float], colors: List[str]) -> None:
        """Plots medians of data as scatter plot"""
        ax.scatter(xs, ys, marker='o', color=colors, s=45, zorder=3, alpha=1)

    @classmethod
    def plot(cls,
             ax,
             data_stack: List[List[float]],
             percentile_ranks: List[float],
             sig_flags: List[bool],
             methods: List[str],
             yaxis_name: Optional[str] = 'EPR (mean)',
             violin_colors: Optional[List[str]] = None,
             x_ticks: Optional[List[str]] = None,
             x_margin: Optional[float] = None,
             legends: Tuple[bool, Optional[List[str]]]=(False,None)
             ) -> None:
        """ Plots violin for the distribution of epr scores
        """
        main_group_n = len(data_stack) / 2  # if it's 16, 8 per each phase
        if violin_colors is None:
            group_repeated_colors = list(np.repeat(cls.group_colors, main_group_n / 2))
            violin_colors = group_repeated_colors + group_repeated_colors
        # - x axis
        if x_ticks is None:
            x_ticks = [i + 1 + (i // main_group_n) for i in
                       range(len(methods))]  # every 8 groups closer to each other
        cls.modify_x_axis(ax, x_ticks, methods, x_margin)
        # - y axis
        cls.modify_y_axis(ax, yaxis_name)
        # - plot medians as scatter plot
        _, medians, _ = np.percentile(data_stack, [25, 50, 75], axis=1)
        cls.plot_median_as_scatter_points(ax, x_ticks, medians, violin_colors)
        # - violin plot
        serif_font()
        bplot = ax.violinplot(np.asarray(data_stack).T, positions=x_ticks, showmeans=True, showextrema=False)

        # - face colors
        for patch_i, patch in enumerate(bplot['bodies']):
            patch.set_facecolor(violin_colors[patch_i])
            patch.set_edgecolor('grey')
            patch.set_alpha(1)
        # - annotate percentile rank
        xs = ax.get_xticks()
        ys = np.max(data_stack, axis=1)
        cls.annotate_percentile_rank(ax, xs, ys, percentile_ranks)
        # - annotate sig sign
        sig_colors = ['red' if item is not None and item < 5 else 'black' for item in percentile_ranks]
        cls.annotate_sig_signs(ax, xs, ys, sig_flags, sig_colors)
        # - plot the legends for the groups
        cls.place_legends_for_groups(ax, legends)
