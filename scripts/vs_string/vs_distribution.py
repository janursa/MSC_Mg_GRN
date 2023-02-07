"""
    Obtains match distribution for two quantiles of .75 and .9
"""
import sys
import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import VS_STRING_DIR, GRN_DIR, ENRICH_DIR, F_DE_data
from scripts.utils import make_title_pretty
from scripts.utils.links import plot_match_counts, determine_sig_signes, compare_network_string, compare_network_string_batch, plot_match_counts_series, create_random_links


if __name__ == '__main__':

    if not os.path.isdir(VS_STRING_DIR):
        os.makedirs(VS_STRING_DIR)

    n_repeat = 10
    methods = ['arbitrary', 'RF', 'ridge', 'portia']
    DE_types = F_DE_data().keys()
    study = 'ctr'

    ncols = 2
    nrows = int(len(DE_types) / ncols)

    links_stack = []
    # - get the links
    for top_quantile in [.75, .9]:
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(4 * ncols, 3 * nrows))
        for idx, DE_type in enumerate(DE_types):
            #- portia
            links_portia = pd.read_csv(os.path.join(GRN_DIR, 'portia', f'links_{DE_type}_{study}.csv'), index_col=False)
            links_stack.append(links_portia)
            mc_portia = compare_network_string(DE_type=DE_type, links=links_portia, top_quantile=top_quantile,
                                   enrich_output_dir=ENRICH_DIR)
            mc_dist_portia = np.repeat(mc_portia, n_repeat)
            # - ridge
            links_ridge = pd.read_csv(os.path.join(GRN_DIR, 'ridge', f'links_{DE_type}_{study}.csv'), index_col=False)
            links_stack.append(links_ridge)
            mc_ridge = compare_network_string(DE_type=DE_type, links=links_ridge, top_quantile=top_quantile,
                                               enrich_output_dir=ENRICH_DIR)
            mc_dist_ridge = np.repeat(mc_ridge, n_repeat)
            #- rf
            links_rf = pd.read_pickle(os.path.join(GRN_DIR, 'RF', DE_type, f'links_pool_{study}.csv'))
            links_stack.append(links_rf)
            mc_dist_rf = compare_network_string_batch(links=links_rf, DE_type=DE_type, top_quantile=top_quantile, enrich_output_dir=ENRICH_DIR)
            # - arbitrary
            links_arbit = create_random_links(links_stack, n=n_repeat)
            mc_dist_arbit = compare_network_string_batch(links=links_arbit, DE_type=DE_type, top_quantile=top_quantile, enrich_output_dir=ENRICH_DIR)

            mc_dist_list = [mc_dist_arbit, mc_dist_rf, mc_dist_ridge, mc_dist_portia]

            sig_signs = determine_sig_signes(mc_dist_list)
            i = int(idx/ncols)
            j = idx%ncols
            ax = axes[i][j]
            ax.set_title(make_title_pretty(DE_type))
            plot_match_counts(ax=ax,data_stack=mc_dist_list, labels=methods, sig_signs=sig_signs)

    #- get the links

        fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_{top_quantile}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_{top_quantile}.pdf'))
        # plt.show()





