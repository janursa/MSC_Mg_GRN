"""
    Obtains top matches for a range of top quantile values.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import ENRICH_DIR, VS_STRING_DIR
from scripts.vs_string.vs_distribution import obtain_links
from scripts.utils.links import plot_match_counts_series, compare_network_string, compare_network_string_batch

if __name__ == '__main__':
    n_repeat = 10
    top_quantiles = np.linspace(.75, .9, 10)
    links_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']
    # - get the links
    links_rf, links_ridge, links_portia, links_random = obtain_links(n_repeat, study='combined')

    ff = lambda links: int(np.mean(compare_network_string_batch(links.copy(), top_quantile, enrich_output_dir=ENRICH_DIR)))
    match_count_list = []
    for top_quantile in top_quantiles:
        match_count_stack = []
        # -- random
        match_count_stack.append(ff(links_random))
        # -- rf
        match_count_stack.append(ff(links_rf))
        # -- ridge
        match_count_stack.append(compare_network_string(links_ridge, top_quantile, enrich_output_dir=ENRICH_DIR))
        # -- portia
        match_count_stack.append(compare_network_string(links_portia, top_quantile, enrich_output_dir=ENRICH_DIR))

        match_count_list.append(match_count_stack)

    match_count_list = np.array(match_count_list).T

    np.savetxt(os.path.join(VS_STRING_DIR,'match_counts_list.csv'), match_count_list, delimiter=",", fmt='%d')

    match_counts_list = np.genfromtxt(os.path.join(VS_STRING_DIR,'match_counts_list.csv'), delimiter=",")
    fig = plot_match_counts_series(match_counts_list,
                                            links_names,
                                            top_quantile_list=top_quantiles)
    fig.savefig(os.path.join(VS_STRING_DIR, 'match_count_trend.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(VS_STRING_DIR, f'match_count_trend.pdf'))
    plt.show()




