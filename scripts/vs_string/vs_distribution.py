"""
    Obtains match distribution for two quantiles of .75 and .9
"""
import sys
import os
import numpy as np
import scipy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import VS_STRING_DIR, GRN_DIR, ENRICH_DIR, F_DE_data
from scripts.utils import make_title_pretty
from scripts.utils.links import compare_network_string, plot_match_counts_series, create_random_links

def determine_sig_signes(datas):
    """ to test sig distribtion from noise links
    Conducts t test to determine whether datas[1:] are significantly different than datas[0], which is ctr
    Datas: Tuple(DataFrame), e.g. [ctr, RF, Ridge, Portia]
    """
    ctr = datas[0] #random
    #- determine p values: compared to ctr
    pvalues = np.array([])
    for data in datas[1:]:
        s, p = scipy.stats.ttest_ind(data, ctr)
        pvalues = np.append(pvalues, p)
    #- determine whether mean distribution is higher than ctr: we only plot higher ones
    increase_flags = np.array((np.mean(datas[1:], axis=1) - np.mean(ctr))>0)
    #- use p values with value flags
    def define_sign(p):
        if p:
            sign = r'$*$'
        else:
            sign=''
        return sign
    flags = (pvalues<0.05)*increase_flags
    sig_signs = ['']+[define_sign(flag) for flag in flags]
    return sig_signs


def plot_match_counts_dist(ax, data_stack, labels, sig_signs):
    matplotlib.rcParams.update({'font.size': 12})

    # fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=(4.7,3.5),
    #     )

    bplot = ax.violinplot(data_stack, showmeans=True, showextrema=False)

    ax.set_ylabel('Number of matched interactions')
    ax.set_xticks(list(range(1,len(labels)+1)))
    ax.set_xticklabels(labels,rotation=0)
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


if __name__ == '__main__':

    if not os.path.isdir(VS_STRING_DIR):
        os.makedirs(VS_STRING_DIR)

    n_repeat = 100
    methods = ['RF', 'ridge', 'portia', 'ensemble', 'arbitrary']
    methods_ordered = ['arbitrary', 'RF', 'ridge', 'portia', 'ensemble']
    DE_types = F_DE_data().keys()
    study = 'ctr'

    ncols = 2
    nrows = int(len(DE_types) / ncols)

    for top_quantile in [.75, .9]: #one plot for each

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(4 * ncols, 3 * nrows))
        for idx, DE_type in enumerate(DE_types):
            links_stack = []
            mc_dist_list = []
            for method in methods:
                if method in ['ensemble', 'portia', 'ridge']:
                    # - single weight
                    links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
                    links_stack.append(links)
                    mc= compare_network_string(DE_type=DE_type, links=links, top_quantile=top_quantile,
                                                       enrich_output_dir=ENRICH_DIR)
                    mc_dist = np.repeat(mc, n_repeat)
                elif method == 'RF':
                    links = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
                    links_stack.append(links)
                    mc_dist = compare_network_string_batch(links=links, DE_type=DE_type, n_repeat=n_repeat, top_quantile=top_quantile, enrich_output_dir=ENRICH_DIR)
                elif method == 'arbitrary':
                    links = create_random_links(links_stack, n=n_repeat)
                    mc_dist = compare_network_string_batch(links=links, DE_type=DE_type, n_repeat=n_repeat,
                                                                 top_quantile=top_quantile,
                                                                 enrich_output_dir=ENRICH_DIR)

                mc_dist_list.append(mc_dist)
            #- change the position of arbit to the first one
            mc_dist_list.insert(0, mc_dist_list[-1])
            del mc_dist_list[-1]


            sig_signs = determine_sig_signes(mc_dist_list)
            i = int(idx/ncols)
            j = idx%ncols
            ax = axes[i][j]
            ax.set_title(make_title_pretty(DE_type))
            plot_match_counts_dist(ax=ax,data_stack=mc_dist_list, labels=methods_ordered, sig_signs=sig_signs)

        fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_{top_quantile}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_{top_quantile}.pdf'))






