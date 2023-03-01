"""
    Obtains top matches for a range of top quantile values and plots different methods vs each other
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR
from scripts.utils.links import plot_match_counts_series, compare_network_string_batch, create_random_links, compare_network_string
from scripts.utils import serif_font, make_title_pretty



if __name__ == '__main__':
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)
    n_repeat = 10
    top_quantiles = np.linspace(.75, .9, 10)

    DE_types = F_DE_data().keys()
    study='ctr'


    ncols = 2
    nrows = int(len(DE_types) / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(4 * ncols, 3 * nrows))


    # - get the links
    for idx, DE_type in enumerate(DE_types):
        if DE_type not in ['late_50']: #TODO: remove
            continue
        mc_stack = []
        methods = []
        for top_quantile in top_quantiles:
            links_stack = []
            # - portia
            links_portia = pd.read_csv(os.path.join(GRN_DIR, 'portia', f'links_{DE_type}_{study}.csv'), index_col=False)
            links_stack.append(links_portia)
            mc_portia = compare_network_string(DE_type=DE_type, links=links_portia, top_quantile=top_quantile,
                                               enrich_output_dir=ENRICH_DIR)
            methods.append('Portia')
            # - ridge
            links_ridge = pd.read_csv(os.path.join(GRN_DIR, 'ridge', f'links_{DE_type}_{study}.csv'), index_col=False)
            links_stack.append(links_ridge)
            mc_ridge = compare_network_string(DE_type=DE_type, links=links_ridge, top_quantile=top_quantile,
                                              enrich_output_dir=ENRICH_DIR)
            methods.append('Ridge')
            # - rf
            links_rf = pd.read_pickle(os.path.join(GRN_DIR, 'RF', f'links_pool_{DE_type}_{study}.csv'))
            links_stack.append(links_rf)
            mc_rf = np.mean(compare_network_string_batch(links=links_rf, DE_type=DE_type, top_quantile=top_quantile, n_repeat=n_repeat,
                                                      enrich_output_dir=ENRICH_DIR))
            methods.append('RF')
            # - arbitrary
            links_arbit, links_arbit_pool = create_random_links(links_stack, n=n_repeat)
            mc_arbit = np.mean(compare_network_string_batch(links=links_arbit_pool, DE_type=DE_type, top_quantile=top_quantile,n_repeat=n_repeat,
                                                         enrich_output_dir=ENRICH_DIR))
            methods.append('Arbitrary')

            mc_list = [mc_portia, mc_ridge, mc_rf, mc_arbit]
            mc_stack.append(mc_list)



        mc_stack = np.array(mc_stack).T


        #--------   plot for all DE_types ---------------------------------

        i = int(idx/ncols)
        j = idx%ncols
        ax = axes[i][j]
        ax.set_title(make_title_pretty(DE_type))
        plot_match_counts_series(mc_stack,
                                 methods,
                                 top_quantile_list=top_quantiles, ax=ax)



    fig.savefig(os.path.join(MODELSELECTION_DIR, f'match_count_trends.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'match_count_trends.pdf'))




