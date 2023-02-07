"""
    Obtains top matches for a range of top quantile values and plots different DE_typs vs each other
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import ENRICH_DIR, VS_STRING_DIR, F_DE_data, GRN_DIR
from scripts.utils.links import plot_match_counts_series, compare_network_string, compare_network_string_batch, create_random_links
from scripts.utils import serif_font, make_title_pretty



if __name__ == '__main__':
    if not os.path.isdir(VS_STRING_DIR):
        os.makedirs(VS_STRING_DIR)
    top_quantiles = np.linspace(.75, .9, 10)

    DE_types = F_DE_data().keys()
    study='ctr'
    methods = ['RF', 'ridge', 'portia']

    ncols = 1
    nrows = 3

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))

    for idx, method in enumerate(methods):
        # - get the links
        mc_stack = [] #for different DE_type
        for DE_type in DE_types:
            mc_list = [] # for different quantiles
            for top_quantile in top_quantiles:
                if method in ['ridge', 'portia']:
                    links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
                    mc = compare_network_string(DE_type=DE_type, links=links, top_quantile=top_quantile,
                                                enrich_output_dir=ENRICH_DIR)
                elif method=='RF':
                    links = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
                    mc = np.mean(compare_network_string_batch(DE_type=DE_type, links=links, top_quantile=top_quantile,
                                                enrich_output_dir=ENRICH_DIR))
                else:
                    raise ValueError('problem')

                mc_list.append(mc)
            mc_stack.append(mc_list)
        mc_stack = np.asarray(mc_stack)

        ax = axes[idx]
        ax.set_title(method)
        plot_match_counts_series(mc_stack,
                                 DE_types,
                                 top_quantile_list=top_quantiles, ax=ax)



    fig.savefig(os.path.join(VS_STRING_DIR, f'mc_trends_method.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(VS_STRING_DIR, f'mc_trends_method.pdf'))




