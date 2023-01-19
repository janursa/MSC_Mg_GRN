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

from imports import OUTPUT_DIR, GRN_DIR
from utils.links import create_random_links, plot_match_counts, choose_top_quantile, compare_network_string, read_write_links, plot_match_counts_series
#- create random links

def compare(links, top_quantile):
    """
    Compare the given links to vs_string
    """
    links_short = choose_top_quantile(links,quantile=top_quantile)
    match_count = compare_network_string(links_short.copy(), OUTPUT_DIR, verbose=False)
    return match_count
def batch_compare(links, top_quantile):
    """
    Compare the given links to vs_string for each weight set in weightpool
    """
    match_counts = []
    weightpool = np.array(links['WeightPool'].values.tolist()).T
    for weight in weightpool:
        links['Weight'] = weight
        match_count = compare(links, top_quantile)
        match_counts.append(match_count)
    return np.array(match_counts)
def determine_sig_signes(datas):
    """
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

def obtain_links(n_repeat):
    # - retreive links for different method
    links_rf = pd.read_pickle(os.path.join(GRN_DIR, 'RF', 'links_pool_ctr.csv'))
    # -- ridge
    links_ridge = read_write_links(study='ctr', method='ridge', mode='read', output_dir=OUTPUT_DIR)
    # -- portia
    links_portia = read_write_links(study='ctr', method='portia', mode='read', output_dir=OUTPUT_DIR)
    # - create random links: run only once and then retreive it
    links_random = create_random_links([links_rf, links_ridge, links_portia], n=n_repeat)
    return links_rf, links_ridge, links_portia, links_random

if __name__ == '__main__':
    n_repeat = 10
    links_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']
    #- get the links
    links_rf, links_ridge, links_portia, links_random = obtain_links(n_repeat)

    for top_quantile in [.75, .9]:
        match_count_stack = []
        # -- random
        match_count_stack.append(batch_compare(links_random.copy(), top_quantile))
        # -- rf
        match_count_stack.append(batch_compare(links_rf.copy(), top_quantile))
        # -- ridge
        match_count_stack.append(np.repeat(compare(links_ridge, top_quantile), n_repeat))
        # -- portia
        match_count_stack.append(np.repeat(compare(links_portia, top_quantile), n_repeat))

        #- print match scores
        for name, item in zip(links_names, match_count_stack):
            print(f'{name}\n \t list size: {item.shape} \n \t mean match score  {np.mean(item)}')

        sig_signs = determine_sig_signes(match_count_stack)
        fig = plot_match_counts(datas=match_count_stack, labels=links_names, sig_signs=sig_signs)

        fig.savefig(os.path.join(OUTPUT_DIR, f'vs_string/match_count_{top_quantile}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(OUTPUT_DIR, f'vs_string/match_count_{top_quantile}.pdf'))
        plt.show()

    # top_quantile_list = np.linspace(.75, .9, num=10)
    #
    # def compare_series():
    #     """
    #     Obtains top matches for a range of top quantile values.
    #     """
    #     match_counts_pool = []
    #     for top_quantile in top_quantile_list:
    #         match_counts = main(links_random, links_rf, links_ridge, links_portia, top_quantile)
    #         match_counts_pool.append(match_counts)
    #     match_counts_list = np.array(match_counts_pool).T
    #     np.savetxt(os.path.join(OUTPUT_DIR,'vs_string','match_counts_list.csv'), match_counts_list, delimiter=",", fmt='%d')
    #
    # compare_series()
    # def plot_series():
    #     match_counts_list = np.genfromtxt(os.path.join(OUTPUT_DIR,'vs_string','match_counts_list.csv'), delimiter=",")
    #     fig = plot_match_counts_series(match_counts_list,
    #                                             links_names,
    #                                             top_quantile_list=top_quantile_list)
    #     fig.savefig(os.path.join(OUTPUT_DIR, f'vs_string/match_count_trend.png'), dpi=300, transparent=True)
    #     fig.savefig(os.path.join(OUTPUT_DIR, f'vs_string/match_count_trend.pdf'))
    #     plt.show()
    # plot_series()



