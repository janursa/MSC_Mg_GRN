"""
description
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *
#- create random links
def create_random_links(links_assembly, n=1000):
    #- TODO: create matrix (using protname)
    links_assembly = [utils.links.nomalize(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j] #flatten
    random_links = links_assembly[0].copy()
    weightpoolvector = []
    for i in range(len(links_assembly[0])):
        weightpoolvector.append(random.sample(weights, n))
    random_links['WeightPool'] = weightpoolvector
    random_links['Weight']= np.mean(weightpoolvector, axis=1)

    random_links.to_pickle(os.path.join(OUTPUT_DIR, 'postprocess',
                                        f'links_random.csv'))
    return random_links
def compare(links, top_quantile):
    """
    Compare the given links to string
    """
    links_short = utils.links.choose_top_quantile(links,quantile=top_quantile)
    match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR, verbose=False)
    return match_count
def batch_compare(links, top_quantile):
    """
    Compare the given links to string for each weight set in weightpool
    """
    match_counts = []
    weightpool = np.array(links['WeightPool'].values.tolist()).T
    for weight in weightpool:
        links['Weight'] = weight
        match_count = compare(links, top_quantile)
        match_counts.append(match_count)
    return match_counts
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
def retrieve_data():
    """
    Read links
    """
    # -- rf
    links_rf = pd.read_pickle(os.path.join(OUTPUT_DIR, 'GRN', 'RF',
                                           f'links_ctr.csv'))  # save with pickle because there is a list of items in data
    # -- ridge
    links_ridge = utils.links.read_write_links('ctr', method='ridge', mode='read_links', OUTPUT_DIR=OUTPUT_DIR)
    # -- portia
    links_portia = utils.links.read_write_links('ctr', method='portia', mode='read_links', OUTPUT_DIR=OUTPUT_DIR)

    return links_rf, links_ridge, links_portia
def main(links_random, links_rf, links_ridge, links_portia, top_quantile=0.75, pool_flag=False, plot_flag=False, links_names=None):
    """
    Main function to run the comparision for all datasets
    """
    # -- ridge
    match_count_ridge = compare(links_ridge, top_quantile)
    # -- portia
    match_count_portia = compare(links_portia, top_quantile)
    if pool_flag:
        match_count_random = batch_compare(links_random.copy(), top_quantile)
        match_count_rf = batch_compare(links_rf.copy(), top_quantile)

        match_count_ridge = [match_count_ridge for i in range(len(match_count_random))]
        match_count_portia = [match_count_portia for i in range(len(match_count_random))]
    else:
        match_count_random = batch_compare(links_random.copy(), top_quantile)
        match_count_rf = batch_compare(links_rf.copy(), top_quantile)
        match_count_random = int(np.mean(match_count_random))
        match_count_rf = int(np.mean(match_count_rf))

        # match_count_rf = np.mean(match_count_random)
    match_counts = match_count_random, match_count_rf, match_count_ridge, match_count_portia
    if plot_flag:
        for label, data in zip(links_names, match_counts):
            print(f'{label}: {np.mean(data)}')
        sig_signs = determine_sig_signes(match_counts)
        fig = utils.links.plot_match_counts(datas=match_counts, labels=links_names, sig_signs=sig_signs)

        fig.savefig(os.path.join(OUTPUT_DIR, f'string/match_count_{top_quantile}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(OUTPUT_DIR, f'string/match_count_{top_quantile}.pdf'))
    return match_counts

if __name__ == '__main__':
    links_names = ['Arbitrary', 'RF', 'Ridge', 'Portia']
    top_quantile_list = np.linspace(.75, .9, num=10)
    links_rf, links_ridge, links_portia = retrieve_data()
    #- create random links: run only once and then retreive it
    links_random = create_random_links([links_rf, links_ridge, links_portia], n=100)
    # links_random = pd.read_pickle(os.path.join(OUTPUT_DIR, 'postprocess', f'links_random.csv'))
    def compare_distribution():
        """
        comapre for top quartile 0.9 and 0.75 with distribution plot
        """
        top_quantile = top_quantile_list[-1]
        main(links_random, links_rf, links_ridge, links_portia, top_quantile, pool_flag=True, plot_flag=True, links_names=links_names)

        top_quantile = top_quantile_list[0]
        main(links_random, links_rf, links_ridge, links_portia, top_quantile, pool_flag=True, plot_flag=True, links_names=links_names)
        plt.show()
    # compare_distribution()
    def compare_series():
        """
        Obtains top matches for a range of top quantile values.
        """
        match_counts_pool = []
        for top_quantile in top_quantile_list:
            match_counts = main(links_random, links_rf, links_ridge, links_portia, top_quantile)
            match_counts_pool.append(match_counts)
        match_counts_list = np.array(match_counts_pool).T
        np.savetxt(os.path.join(OUTPUT_DIR,'string','match_counts_list.csv'), match_counts_list, delimiter=",", fmt='%d')

    # compare_series()
    def plot_series():
        match_counts_list = np.genfromtxt(os.path.join(OUTPUT_DIR,'string','match_counts_list.csv'), delimiter=",")
        fig = utils.links.plot_match_counts_series(match_counts_list,
                                                links_names,
                                                top_quantile_list=top_quantile_list)
        fig.savefig(os.path.join(OUTPUT_DIR, f'string/match_count_trend.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(OUTPUT_DIR, f'string/match_count_trend.pdf'))
        plt.show()
    plot_series()



