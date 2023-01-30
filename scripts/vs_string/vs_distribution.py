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

from scripts.imports import VS_STRING_DIR, GRN_DIR, ENRICH_DIR
from scripts.utils.links import create_random_links, plot_match_counts, determine_sig_signes, compare_network_string, compare_network_string_batch, read_write_links, plot_match_counts_series


def obtain_links(n_repeat, study):
    # - retreive links for different method
    links_rf = pd.read_pickle(os.path.join(GRN_DIR, 'RF', 'links_pool_ctr.csv'))
    # -- ridge
    links_ridge = read_write_links(study=study, method='ridge', mode='read', output_dir=GRN_DIR)
    # -- portia
    links_portia = read_write_links(study=study, method='portia', mode='read', output_dir=GRN_DIR)
    # - create random links: run only once and then retreive it
    links_random = create_random_links([links_rf, links_ridge, links_portia], n=n_repeat)
    return links_rf, links_ridge, links_portia, links_random

if __name__ == '__main__':

    if not os.path.isdir(VS_STRING_DIR):
        os.makedirs(VS_STRING_DIR)
    n_repeat = 10
    links_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']
    #- get the links
    links_rf, links_ridge, links_portia, links_random = obtain_links(n_repeat, study='ctr')

    for top_quantile in [.75, .9]:
        match_count_stack = []
        # -- random
        match_count_stack.append(compare_network_string_batch(links_random.copy(), top_quantile, enrich_output_dir=ENRICH_DIR))
        # -- rf
        match_count_stack.append(compare_network_string_batch(links_rf.copy(), top_quantile, enrich_output_dir=ENRICH_DIR))
        # -- ridge
        match_count_stack.append(np.repeat(compare_network_string(links_ridge, top_quantile, enrich_output_dir=ENRICH_DIR), n_repeat))
        # -- portia
        match_count_stack.append(np.repeat(compare_network_string(links_portia, top_quantile, enrich_output_dir=ENRICH_DIR), n_repeat))

        #- print match scores
        # for name, item in zip(links_names, match_count_stack):
        #     print(f'{name}\n \t list size: {item.shape} \n \t mean match score  {np.mean(item)}')

        sig_signs = determine_sig_signes(match_count_stack)
        fig = plot_match_counts(datas=match_count_stack, labels=links_names, sig_signs=sig_signs)

        fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_{top_quantile}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_{top_quantile}.pdf'))
        plt.show()





